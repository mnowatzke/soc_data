library(dplyr)

# Author: Jarad Niemi
# Date:   2024-01-21
# Purpose: Construct covariance matrices
#
# Details: The purpose of this script is to describe how to construct
#   covariance matrices in order to utilize these constructed covariance
#   matrices in Bayesian estimation procedures, e.g. DREAM.
#
#   The script will attempt to reconstruct covariance matrices used in
#   Mathers et. al. (2023) equation (1) explained in words in the second
#   paragraph of Page 6.
#
#   The script will then attempt to provide code that can be used by SEC to
#   incorporate appropriate observation dependencies.
#
#   All positive parameters that may be estimated are included as their logarithms
#   so that optimization or MCMC can be performed on an unrestricted scale.



#' Construct independent covariance
#'
#' Construct diagonal covariance matrix from logarithm of the variance and
#' number of observations.
#'
#' @param log_sigma2 numeric, the logarithm of the variance
#' @param n integer, the number of observations
#' @return an n x n, (diagonal) covariance matrix
#' @examples
#' construct_independent_covariance(0.1)
#'
construct_independent_covariance <- function(log_sigma2, n = 10) {
  exp(log_sigma2)*diag(n)
}

#' Temporal independent covariance
#'
#' In Mathers et. al. (2023), the used an independent (residual) covariance
#' matrix that dependended on t: number of years since the first SOC measurement.
#' This function will create that type of covariance matrix.
#'
#' @param log_sigma2 numeric, logarithm of sigma2
#' @param log_v numeric, logarithm of nu
#' @param t numeric (integer), number of years since first SOC measurement
#' @return length(t) x length(t) covariance matrix
#' @examples
#' construct_temporal_independent_covariance(log(1), log(0.001), 1:10)
#'
construct_temporal_independent_covariance <- function(log_sigma2, log_v, t) {
  diag(exp(log_sigma2) * exp(2*t*exp(log_v)))
}


#' Construct covariance for random effect matrix
#'
#' In Mathers et. al. (2023) each site had a random effect and observations from
#' the same site are correlated.
#'
#' @param log_sigma2_site numeric, logarithm of the variance of the site
#' @param site factor (or character), indicating siteIDs
#' @return length(site) x length(site) covariance matrix
#' @examples
#' # site random effect: 3 sites, each measured twice
#' construct_random_effect_covariance(log(1), rep(c("A","B","C"), each = 2))
#'
#' # year-within-site random effect: 2 sites, each measured twice in two different years
#' d <- expand.grid(rep = 1:2, year = 2021:2022, site = c("A","B"))
#' d$siteyear <- paste0(d$site, d$year)
#' construct_random_effect_covariance(log(1), d$siteyear)
#'
construct_random_effect_covariance <- function(log_sigma2, random_effect) {
  d <- data.frame(random_effect = factor(random_effect))
  m <- as.matrix(model.matrix(~0+random_effect, data = d))
  exp(log_sigma2) * m %*% t(m)
}


# Putting these together, I believe the covariance structure for Mathers et. al.
# (2023) can be constructed using this function
construct_mathers_covariance <- function(params, d) {
  # Site random effect
  construct_random_effect_covariance(
    params[1],
    d$site) +

    # Year-within-site random effect
    construct_random_effect_covariance(
      params[2],
      paste0(
        d$site,
        d$year)) +

    # Time-based residual variance
    construct_temporal_independent_covariance(
      params[3],
      params[4],
      d$years_since_SOC_measurement
    )
}


set_cov_cols <- function(studies_df) {
  studies_df$site_id <- as.numeric(factor(studies_df$location))
  studies_df$yrs_btwn_meas <- studies_df$trt_yr_end - studies_df$trt_yr_st
  studies_df$rep_meas <- 0
  for (id in unique(studies_df$site_id)) {
    sub_df <- subset(studies_df, site_id == id)
    n_meas <- length(unique(sub_df$trt_yr_end))
    idx <- which(studies_df$site_id == id)
    studies_df$rep_meas[idx] <- n_meas
  }
  return(studies_df)
}



constr_mather_cov <- function(params, site_ids, obs_years, yrs_btwn_meas) {
  # Site random effect
  construct_random_effect_covariance(
    params[1],
    site_ids) +

    # Year-within-site random effect
    construct_random_effect_covariance(
      params[2],
      paste0(
        site_ids,
        obs_years)) +

    # Time-based residual variance
    construct_temporal_independent_covariance(
      params[3],
      params[4],
      yrs_btwn_meas
    )
}