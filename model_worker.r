#### MCMC sampler and predictions for APSIM model
#### Uncertainty estimates for changes in SOC and yield
library(optparse)
library(BayesianTools)
library(configr)
library(xml2)
library(tidyr)
library(dplyr)
library(reticulate)
library(RPostgreSQL)
library(ggplot2)
library(yaml)
library(uuid)
library(stringr)
library(parallel)
source("src/parse_out.r")
source("src/edit_param.r")
source("src/create_plots.r")
source("src/get_config.r")
source("src/get_sim_delta.r")

### Command line args
option_list = list(
  make_option(c("-t", "--trainfolder"), type="character", default=NULL, 
      help="folder with training configuration files", metavar="character"),
  make_option(c("-v", "--valfolder"), type="character", default=NULL, 
      help="folder with validation configuration files", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# set cores to use
# tot_cores <- detectCores()
# n_cores <- as.integer(tot_cores - 2)

#set our configuration training and testing folders as variables
train_config_dir <- opt$trainfolder
val_config_dir <- opt$valfolder
#read in the bayes model settings
model_config <- read.config('src/model_config.toml')
model_params <- model_config$model$calib_params

# set python environment and run_apsim script
# use_python("C:/vmd53/env/Scripts/python.exe")
use_condaenv("sand")
source_python('src/run_apsim.py')

# get path to our vmd53 folder
vmd53_folder <- getwd()
# vmd53_folder <- dirname(getwd())

### get our params of interest by rainfed and irrigated regions
params_df <- read.csv(file.path(vmd53_folder, 'parameters', 'calib_params.csv'))
params_df <- dplyr::filter(params_df, param %in% model_params)
message('Using parameters:')
message(cat(model_params))

obs_studies <- read.csv(file.path(vmd53_folder, 'data', 'obs', 'tillage_papers.csv'))
# if observation is in Mg ha-1, convert to kg ha-1
obs_studies$value[obs_studies$unit == 'Mg ha-1'] <- obs_studies$value * 1000
obs_studies$yrly_delta[obs_studies$unit == 'Mg ha-1'] <- obs_studies$yrly_delta * 1000
obs_studies$unit[obs_studies$unit == 'Mg ha-1'] <- 'kg ha-1'
setDT(obs_studies)

if (!is.null(train_config_dir)) {
    #start with training and get all config files
    train_config_files <- list.files(train_config_dir, pattern="\\.toml$", full.names = TRUE)
    if (length(train_config_files) == 0) {
        stop(paste0("No configuration files found in ", train_config_dir), call.=FALSE)
    }
    message(paste0('Found ', length(train_config_files), ' configuration files for training.'))
    train_config_df <- get_config_df(train_config_files)
    #check to see if commodity, rotation, and mgmt intervention are the same
    ifelse (
        length(unique(train_config_df[,'study_commodity'])) == 1,
        commodity <- train_config_df[,'study_commodity'][1],
        stop('Commodity is not the same for all training config files')
    )
    ifelse (
        length(unique(train_config_df[,'study_rot'])) == 1,
        rotation <- train_config_df[,'study_rot'][1],
        stop('Commodity is not the same for all training config files')
    )
    ifelse (
        length(unique(train_config_df[,'study_interv'])) == 1,
        mgmt_interv <- train_config_df[,'study_interv'][1],
        stop('Intervention is not the same for all training config files')
    )
    message('Training configuration files match for commodity ', commodity, ', rotation ', rotation, ', and intervention ', mgmt_interv)

    # set our xml and met folders
    met_folder <- file.path(vmd53_folder, 'data', 'apsim_files', 'train', mgmt_interv, rotation, 'met_files')
    xml_folder <- file.path(vmd53_folder, 'data', 'apsim_files', 'train', mgmt_interv, rotation, 'crop_xml')

    apsim_files_dir <- file.path(vmd53_folder, 'data', 'apsim_files', 'train', mgmt_interv, rotation)
    #copy our apsim files to run
    #check if there are any apsim files with leading _ as these are the copied
    copied_apsim_files <- list.files(apsim_files_dir, pattern="^_[A-Za-z]+_.*\\.apsim$")
    original_apsim_files <- list.files(apsim_files_dir, pattern="^[A-Za-z]+_.*\\.apsim$")
    if (length(original_apsim_files) == 0) {
        stop(cat("Missing apsim files for directory ", apsim_files_dir))
    }
    # if copied sims already exist, overwrite them
    if (length(copied_apsim_files) != 0) {
        copied_apsim_paths <- list.files(apsim_files_dir, pattern="^_[A-Za-z]+_.*\\.apsim$", full.names = TRUE)
        lapply(original_apsim_files, function(x) file.copy(file.path(apsim_files_dir, x), file.path(apsim_files_dir, paste0('_', x)), overwrite=TRUE))
    } else {
        #else copy them 
        lapply(original_apsim_files, function(x) file.copy(file.path(apsim_files_dir, x), file.path(apsim_files_dir, paste0('_', x)), overwrite=TRUE))
        copied_apsim_paths <- list.files(apsim_files_dir, pattern="^_[0-9]+_.*\\.apsim$", full.names = TRUE)
        # copied_apsim_paths <- file.path(apsim_files_dir, copied_apsim_files)
    }
    n_apsim_files = length(original_apsim_files)
    message(paste0("Using runs for ", n_apsim_files, " apsim files/soil mukeys."))
        # set cores to use
    tot_cores <- detectCores()
    if (n_apsim_files >= tot_cores) {
        n_cores <- as.integer(tot_cores - 2)
    } else {
        n_cores <- n_apsim_files
    }
    # get our commodity yield var of interest by matching with the commodity in config
    crop_yield_df = data.frame(crops=c('soybean', 'maize', 'wheat', 'rye'),
                                yld_vars=c('soy_yield', 'maize_yield', 'wheat_yield', 'rye_biomass'))
    yld_comm = crop_yield_df$yld_vars[match(commodity, crop_yield_df$crops)]

    # run all of our simulations once to see how it performs
    run_all_sims = try(run_apsim_main(apsim_files_dir, n_cores=n_cores))
    if(inherits(run_all_sims, "try-error")) {
            print('Error running simulation.')
    }
    # create dataframe so save pre-optimized oc delta
    pre_optim_df <- data.frame(
        study=character(),
        site=character(),
        rot=character(),
        comm=character(),
        interv=character(),
        yr_st=numeric(),
        yr_end=numeric(),
        lyr_top=numeric(),
        lyr_btm=numeric(),
        sim=numeric(),
        obs=numeric()
    )

    ll_df <- data.frame(
        study=character(),
        site=character(),
        rot=character(),
        comm=character(),
        interv=character(),
        yr_st=numeric(),
        yr_end=numeric(),
        lyr_top=numeric(),
        lyr_btm=numeric(),
        sim=numeric(),
        obs=numeric()
    )
    # loop through each of the config files
    for (row in 1:nrow(train_config_df)) {
        study_name = train_config_df[row, 'study_name']
        study_site = train_config_df[row, 'study_site']
        study_rot = train_config_df[row, 'study_rot']
        study_commodity = train_config_df[row, 'study_commodity']
        study_interv = train_config_df[row, 'study_interv']
        study_yr_st = train_config_df[row, 'study_yr_st']
        study_yr_end = train_config_df[row, 'study_yr_end']
        # get individual study info
        study_obs = obs_studies[
            study == study_name &
            location == study_site &
            rotation == study_rot &
            treatment == study_interv &
            commodity == study_commodity #&
            # trt_yr_st == study_yr_st &
            # trt_yr_end == study_yr_end 
        ]
        if (length(study_obs) == 0) {
            stop(cat('Could not find matching observations for study ', study_name, ', location ', study_site,
            ', rotation ', study_rot, ', treatment ', study_interv, ', and commodity', study_commodity))
        }
        # read in out files
        study_out_pattern = paste0(
            "study_", study_name, '_site_', study_site, '_interv_',
            study_interv, '_rot_', study_rot, '_[[:alnum:]_-]+\\.out$'
        )
        out_files <- list.files(path=apsim_files_dir, pattern=study_out_pattern, full.names=TRUE)
        for (obs_row in 1:nrow(study_obs)) {
            oc_top = as.integer(study_obs[obs_row, 'lyr_top'])
            oc_btm = as.integer(study_obs[obs_row, 'lyr_btm'])
            study_yr_st = as.integer(study_obs[obs_row, 'trt_yr_st'])
            study_yr_end = as.integer(study_obs[obs_row, 'trt_yr_end'])
            obs_delta = as.numeric(study_obs[obs_row, 'yrly_delta'])
            # obs_delta = obs_delta$yrly_delta
            pre_sum_df = sum_multi_sims(
                out_files,
                study_yr_st,
                study_yr_end,
                oc_top,
                oc_btm
            )
            pre_hb_delta = round(mean(pre_sum_df$hb_delta_yr), 2)
            if(is.na(pre_hb_delta)) print(paste0('OC delta not valid for ', study_obs$study, ' ',
                                                    study_obs$location))
            if(is.numeric(pre_hb_delta) != TRUE) {
                pre_hb_delta = as.numeric(pre_hb_delta)
            }
            pre_optim_df[nrow(pre_optim_df) + 1,] = list(
                study_name, study_site, study_rot, study_commodity,
                study_interv, study_yr_st, study_yr_end, oc_top, oc_btm, pre_hb_delta, obs_delta
            )
            ll_df[nrow(ll_df) + 1,] = list(
                study_name, study_site, study_rot, study_commodity,
                study_interv, study_yr_st, study_yr_end, oc_top, oc_btm, pre_hb_delta, obs_delta
            )
        }
        # TODO check if study config crops match .apsim yield output
        # qc_yld_out(study_config, pre_crop_ylds)
    }
    pre_err <- pre_optim_df$obs - pre_optim_df$sim 
    message(cat("Pre-optimized error is ", pre_err))
    # note: running for 7 cc cp fields the sd was 428.60 kg ha-1
    pre_sd <- sd(pre_err)
    message(cat("Pre-optimized standard deviation is ", pre_sd))
    model_out_dir <- file.path(vmd53_folder, "data", "model_out", mgmt_interv, rotation)
    if (!dir.exists(model_out_dir)) {
        message(paste0('Directory does not exist for ', model_out_dir, ', creating directory.'))
        dir.create(model_out_dir, recursive = TRUE)
    }
    save.image(file=file.path(model_out_dir, paste0('pre_optim_', mgmt_interv, '_', rotation, '_sd.RData')))
    # get the .sim file paths for editing arameters in the .sim files
    # easiest to get these filenames after running the apsim files once already
    sim_paths <- list.files(apsim_files_dir, pattern="^study.*\\.sim$", full.names = TRUE)
    
    n_params <- nrow(params_df)
    start_matrix <- matrix(ncol=n_params, nrow=3)
    colnames(start_matrix) <- params_df$param
    for (i in 1:nrow(params_df)) {
        lower = params_df[i, 'lower']
        upper = params_df[i, 'upper']
        middle = (lower + upper)/2
        start_vals = c(lower, upper, middle)
        start_matrix[,i] <- start_vals
    }
    loglik_apsim <- function(param_vals) {
        # set parameters in dataframe to input values
        params_df$param_vals <- param_vals
        # read in apsim file and write with new param values for params in .apsim file
        remove_old_files(apsim_files_dir)
        apsim_params <- dplyr::filter(params_df, param_src %in% c('apsim'))
        if (nrow(apsim_params) != 0) {
            for (apsim_path in copied_apsim_paths) {
                apsim_xml <- read_xml(apsim_path)
                set_params(apsim_params, sim_xml=apsim_xml)
                write_xml(apsim_xml, apsim_path)
            }
        }
        # convert .apsim to .sim files
        convert_all_apsim = try(convert_all_apsim_to_sim(copied_apsim_paths, num_cores=n_cores))
        if(inherits(convert_all_apsim, "try-error")) {
                print('Error converting .apsim files to .sim files.')
                return(NA)
        }
        # edit our .sim params
        sim_params <- dplyr::filter(params_df, param_src %in% c('sim'))
        if (nrow(sim_params) != 0) {
            for (sim_path in sim_paths) {
                sim_xml <- read_xml(sim_path)
                set_params(sim_params, sim_xml=sim_xml)
                write_xml(sim_xml, sim_path)
            }
        }
        # run all .sim files
        run_all_sims = try(run_many_sims(sim_paths, num_cores=n_cores))
        if(inherits(run_all_sims, "try-error")) {
                print('Error running simulation.')
                return(NA)
        }
        remove_tmp_files(apsim_files_dir)
        ll_df <- get_sim_diff_delta(ll_df, train_config_df, obs_studies, apsim_files_dir)
        err <- ll_df$obs - ll_df$sim
        #calculate the likelihood
        #TODO 
        # get new parameters for covariance matrix into param_vals
        # construct the cov matrix
        # NOTE if we run into errors with cov matrix not being > 0,
        # (error not pivitive definite etc), add diagonal matrix 0.00001,
        # possibly can use machine precision sqrt(.Machine$double.eps) for this
        # construct new priors for additional cov terms/params
        # mvtnrom::dmvnorm(ll_df$obs, mean=ll_df$sim, sigma=cov_matrix, log=TRUE)
        loglik <- stats::dnorm(err, mean = 0, sd = pre_sd, log = TRUE)
        return(sum(loglik))
    }
    prior <- createTruncatedNormalPrior(mean = params_df$default,
                                    sd = params_df$sd,
                                    lower = params_df$lower,
                                    upper = params_df$upper)

    bayes_setup <- BayesianTools::createBayesianSetup(loglik_apsim, 
                                                        prior,
                                                        names = params_df$param)

    # sensitivity analysis
    # plotSensitivity(bayes_setup)
    # #TODO check for convergence
    # coda::gelman.diag(chain)
    mcmc_iterations <- model_config$model$iterations
    mcmc_nchains <- model_config$model$n_chains
    mcmc_burnin <- model_config$model$burnin
    mcmc_sampler <- model_config$model$sampler
    settings <- list(iterations = mcmc_iterations, nrChains = mcmc_nchains, startValue=start_matrix)
    message(cat(paste0('Running MCMC for ', mcmc_iterations, ' iterations and ', mcmc_nchains, ' chains.')))
    # markovs and chains
    mcmc_out <- runMCMC(bayesianSetup = bayes_setup, sampler = mcmc_sampler, settings = settings)
    save.image(file=file.path(model_out_dir, paste0('train_', mgmt_interv, '_', rotation, '_mcmcout.RData')))
    save(mcmc_out, file=file.path(model_out_dir, paste0('train_', mgmt_interv, '_', rotation, '_mcmcout.rds')))
    # plot images and save to dir
    img_out_dir <- file.path(model_out_dir, "imgs")
    if (!dir.exists(img_out_dir)) {
        dir.create(img_out_dir)
    }
    message('Saving MCMC plots for mcmc out, parameter marginals, and covariance matrix.')
    pdf(file=file.path(img_out_dir, 'mcmc_plots.pdf'))
    plot(mcmc_out)
    marginalPlot(mcmc_out)
    if (nrow(params_df) > 1) {
        correlationPlot(mcmc_out)
    }
    dev.off()
    #NOTE if ggsave fails may need to reinstall package 'ragg'
    # get diagnostic info on mcmc
    gel_diag <- gelmanDiagnostics(mcmc_out, plot=T)
    mcmc_sum <- summary(mcmc_out)
    message('Summary for MCMC')
    print(mcmc_sum)
    message('Gelman diagnostics for MCMC sampler.')
    print(gel_diag)
}

if (!is.null(val_config_dir)) {
    # get configuration toml files
    val_config_files <- list.files(val_config_dir, pattern="\\.toml$", full.names = TRUE)
    if (length(val_config_files) == 0) {
        stop(paste0("No validation configuration files found in ", val_config_dir), call.=FALSE)
    }
    message(paste0('Found ', length(val_config_files), ' configuration files for validationg.'))
    # check to see if configuration files have same commodity, rotation, and intervention
    val_config_df <- get_config_df(val_config_files)
    ifelse (
        length(unique(val_config_df[,'study_commodity'])) == 1,
        commodity <- val_config_df[,'study_commodity'][1],
        stop('Commodity is not the same for all validation config files')
    )
    ifelse (
        length(unique(val_config_df[,'study_rot'])) == 1,
        rotation <- val_config_df[,'study_rot'][1],
        stop('Commodity is not the same for all validation config files')
    )
    ifelse (
        length(unique(val_config_df[,'study_interv'])) == 1,
        mgmt_interv <- val_config_df[,'study_interv'][1],
        stop('Intervention is not the same for all validation config files')
    )
    message('Validation configuration files match for commodity ', commodity, ', rotation ', rotation, ', and intervention ', mgmt_interv)
    # if just validating, read in the Rdata file.
    model_out_dir <- file.path(vmd53_folder, "data", "model_out", mgmt_interv, rotation)
    if (!dir.exists(model_out_dir)) {
        message(paste0('Creating directory at ', model_out_dir))
        dir.create(model_out_dir, recursive = TRUE)
    }
    mcmc_out_file = file.path(model_out_dir, paste0('train_', mgmt_interv, '_', rotation, '_mcmcout.rds'))
    # read in Rdata file if exists
    if (file.exists(mcmc_out_file)) {
        load(file=mcmc_out_file)
    }

    # get apsim files to run for validation and copy them if needed
    val_files_dir <- file.path(vmd53_folder, 'data', 'apsim_files', 'val', mgmt_interv, rotation)
    vcopied_apsim_files <- list.files(val_files_dir, pattern="^_[A-Za-z]+_.*\\.apsim$")
    voriginal_apsim_files <- list.files(val_files_dir, pattern="^[A-Za-z]+_.*\\.apsim$")
    if (length(voriginal_apsim_files) == 0) {
        stop(cat("Missing apsim files for directory ", val_files_dir))
    }
    # if copied sims already exist, overwrite them
    if (length(vcopied_apsim_files) != 0) {
        vcopied_apsim_paths <- list.files(val_files_dir, pattern="^_[0-9]+_.*\\.apsim$", full.names = TRUE)
        lapply(voriginal_apsim_files, function(x) file.copy(file.path(val_files_dir, x), file.path(val_files_dir, paste0('_', x)), overwrite=TRUE))
    } else {
        #else copy them 
        lapply(voriginal_apsim_files, function(x) file.copy(file.path(val_files_dir, x), file.path(val_files_dir, paste0('_', x)), overwrite=TRUE))
        vcopied_apsim_paths <- list.files(val_files_dir, pattern="^_[0-9]+_.*\\.apsim$", full.names = TRUE)
        # vcopied_apsim_paths <- file.path(val_files_dir, vcopied_apsim_files)
    }
    n_apsim_files = length(voriginal_apsim_files)
    message(paste0("Using runs for ", n_apsim_files, " apsim files/soil mukeys."))
    # set cores to use
    tot_cores <- detectCores()
    if (n_apsim_files >= tot_cores) {
        n_cores <- as.integer(tot_cores - 2)
    } else {
        n_cores <- n_apsim_files
    }
    # get our commodity yield var of interest by matching with the commodity in config
    crop_yield_df <- data.frame(crops=c('soybean', 'maize', 'wheat', 'rye'),
                                yld_vars=c('soy_yield', 'maize_yield', 'wheat_yield', 'rye_biomass'))
    yld_comm <- crop_yield_df$yld_vars[match(commodity, crop_yield_df$crops)]

    # run all of our simulations once to see how it performs
    run_val_sims = try(run_apsim_main(val_files_dir, n_cores=n_cores))
    if(inherits(run_val_sims, "try-error")) {
            print('Error running simulations.')
    }
    sim_paths <- list.files(val_files_dir, pattern="^study.*\\.sim$", full.names = TRUE)
    post_pred_df <- data.frame(
        study=character(),
        site=character(),
        rot=character(),
        comm=character(),
        yr_st=integer(),
        yr_end=integer(),
        interv=character(),
        sim=numeric(),
        obs=numeric()
    )
    for (row in 1:nrow(params_df)) {
        optim_param = params_df[row, 'param']
        post_pred_df[,optim_param] <- numeric()
    }
    post_pred_df[,'avg_pred'] <- numeric()
    # post_pred_df[,'err'] <- numeric()
    # post_pred_df[,'mean_err'] <- numeric()
    create_predictions <- function(param_vals) {
        # set parameters in dataframe to input values
        params_df$param_vals <- param_vals
        # read in apsim file and write with new param values for params in .apsim file
        remove_old_files(val_files_dir)
        apsim_params <- dplyr::filter(params_df, param_src %in% c('apsim'))
        if (nrow(apsim_params) != 0) {
            for (apsim_path in vcopied_apsim_paths) {
                apsim_xml <- read_xml(apsim_path)
                set_params(apsim_params, sim_xml=apsim_xml)
                write_xml(apsim_xml, apsim_path)
            }
        }
        # convert .apsim to .sim files
        convert_all_apsim = try(convert_all_apsim_to_sim(vcopied_apsim_paths, num_cores=n_cores))
        if(inherits(convert_all_apsim, "try-error")) {
                print('Error converting .apsim files to .sim files.')
        }
        # edit our .sim params
        sim_params <- dplyr::filter(params_df, param_src %in% c('sim'))
        if (nrow(sim_params) != 0) {
            for (sim_path in sim_paths) {
                sim_xml <- read_xml(sim_path)
                set_params(sim_params, sim_xml=sim_xml)
                write_xml(sim_xml, sim_path)
            }
        }
        # run all .sim files
        run_all_sims = try(run_many_sims(sim_paths, num_cores=n_cores))
        if(inherits(run_all_sims, "try-error")) {
                print('Error running simulation.')
                return(NA)
        }
        remove_tmp_files(val_files_dir)

        pred_df = data.frame(
            study=character(),
            site=character(),
            rot=character(),
            comm=character(),
            yr_st=integer(),
            yr_end=integer(),
            lyr_top=numeric(),
            lyr_btm=numeric(),
            interv=character(),
            sim=numeric(),
            obs=numeric()
        )
        for (row in 1:nrow(params_df)) {
            optim_param = params_df[row, 'param']
            pred_df[,optim_param] <- numeric()
        }
        # do some summary statistics with .out files

        for (row in 1:nrow(val_config_df)) {
            # get individual study info
            study_name = val_config_df[row, 'study_name']
            study_site = val_config_df[row, 'study_site']
            study_rot = val_config_df[row, 'study_rot']
            study_commodity = val_config_df[row, 'study_commodity']
            study_interv = val_config_df[row, 'study_interv']
            study_yr_st = val_config_df[row, 'study_yr_st']
            study_yr_end = val_config_df[row, 'study_yr_end']
            # get our observations
            study_obs = obs_studies[
                study == study_name &
                location == study_site &
                rotation == study_rot &
                treatment == study_interv &
                commodity == study_commodity &
                trt_yr_st == study_yr_st &
                trt_yr_end == study_yr_end 
            ]
            if (nrow(study_obs) == 0) {
                next
                cat("No observations found for study ", study_name, " ", study_site, " ", study_rot, " ", study_interv, " ", study_commodity)
            }
            # do some summary statistics with .out files
            # read in out files
            study_out_pattern = paste0(
                "study_", study_name, '_site_', study_site, '_interv_',
                study_interv, '_rot_', study_rot, '_[[:alnum:]_-]+\\.out$'
            )
            out_files = list.files(path=val_files_dir, pattern=study_out_pattern, full.names=TRUE)
            for (obs_row in 1:nrow(study_obs)) {
                oc_top = as.integer(study_obs[obs_row, 'lyr_top'])
                oc_btm = as.integer(study_obs[obs_row, 'lyr_btm'])
                obs_delta = round(study_obs[obs_row, 'yrly_delta'], 2)
                obs_delta = obs_delta$yrly_delta
                sum_df = sum_multi_sims(
                    out_files,
                    study_yr_st,
                    study_yr_end,
                    oc_top,
                    oc_btm
                )
                hb_delta = round(mean(sum_df$hb_delta_yr), 2)
                if(is.na(hb_delta)) print(paste0('OC delta not valid for ', study_obs$study, ' ',
                                                    study_obs$location))
                pred_df_list = list(study_name, study_site,
                study_rot, study_commodity, study_yr_st, study_yr_end, oc_top, oc_btm,
                study_interv, hb_delta, obs_delta
                )
                for (row in 1:nrow(params_df)) {
                    param_val = params_df[row, 'param_vals']
                    pred_df_list <- append(pred_df_list, param_val)
                }
                # append to pred df
                pred_df[nrow(pred_df) + 1,] = pred_df_list
            }
        }
        predic = round(mean(pred_df$sim), 2)
        pred_df$avg_pred = predic
        # pred_df$err = pred_df$obs - pred_df$sim
        # pred_df$mean_err = round(mean(err), 2)
        post_pred_df <<- rbind(post_pred_df, pred_df)
        return(predic)
    }
    img_out_dir <- file.path(model_out_dir, "imgs")
    message(paste0('Sampling parameters from MCMC out for new predictions. See plot at ', file.path(img_out_dir, "post_pred_params.png")))
    #TODO set this in config or cmd
    num_samples <- model_config$pred$iterations
    param_matrix <- getSample(mcmc_out, numSamples = num_samples)
    # check for number of samples to draw, 1k or n, whichever is less
    num_samples <- min(1000, nrow(param_matrix))
    post_param_df <- data.frame(param_matrix)
    pdf(file=file.path(img_out_dir, 'post_param_samples.pdf'))
    create_dense_plots(post_param_df, title='Distributions of sampled parameters')
    dev.off()

    #### make new predictions from MCMC params ####
    message(paste0('Getting post param predictions for ', num_samples, ' samples.'))
    pred_dist <- getPredictiveDistribution(param_matrix, model=create_predictions, numSamples = num_samples)
    save(pred_dist, file=file.path(model_out_dir, paste0(mgmt_interv, '_', rotation, '_', commodity, '_pred_dist.rds')))
    save(post_pred_df, file=file.path(model_out_dir, paste0(mgmt_interv, '_', rotation, '_', commodity, '_post_pred.rds')))
    message('Finished with predictions.')
    pred_dist_df <- data.frame(pred_dist)
    pred_dist_df <- cbind(param_matrix, pred_dist_df)
    med_oc_delta <- median(pred_dist_df$med_oc_delta)
    # med_yield <- median(pred_dist_df$yield)
    # message(paste0('Median post predicted yield is ', med_yield))
    message(paste0('Median post oc delta is ', med_oc_delta))
    # save as plot
    pdf(file=file.path(img_out_dir, 'post_pred.pdf'))
    create_dense_plots(pred_dist_df, title='Distribution of posterior OC change predictions')
    dev.off()
    # get our 95% credible intervals
    cred_intervs <- get_cred_int(pred_dist, quantiles=c(0.05, 0.95))
    message(paste0('SOC delta 90% credible interval: ', cred_intervs[1,1], ' to ', cred_intervs[2,1], ' kg ha-1'))
    message(paste0('Saving all output to ', model_out_dir, ' for ', commodity, ' ', rotation, ' ', mgmt_interv))
    save.image(file=file.path(model_out_dir, paste0(mgmt_interv, '_', rotation, '_', commodity, '_preds.RData')))
}