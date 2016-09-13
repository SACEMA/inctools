# Created by and Copyright (C) 2015-2016 Stefano Ongarello (FIND). 
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

#' Calculates the required number of subjects in a survey in case 1
#'
#' @param TIME in days
#' @param MDRI in days
#' @param mdrihatcov covariance associated with MDRI estimate
#' @param frrhat FRR estimate
#' @param frrhatcov covariance associated with FRR estimate
#' @param inc_1 incidence in survey 1
#' @param inc_2 incidence in survey 2
#' @param p_pos_1 prevalence in survey 1
#' @param p_pos_2 prevalence in survey 2
#' @param power required power
#' @param alpha level of significance
#' @param DE_prev_1 DE prevalence 1
#' @param DE_prev_2 DE prevalence 2
#' @param DE_RgivenTested_1 DE prevalence given tested 1
#' @param DE_RgivenTested_2 DE prevalence given tested 2
#' @return sample_size the corresponding sample siye as required by inputs
#' @param rec_test_coverage_1 recent test coverage 1
#' @param rec_test_coverage_2 recent test coverage 2
ss_calc_case1 <- function(TIME = 730, MDRI = 200,
                             mdrihatcov = 0.05, frrhat = 0.01, frrhatcov = 0.20,
                             inc_1 = 0.05, inc_2 = 0.03, p_pos_1 = 0.20, p_pos_2 = 0.15,
                             power = 0.8,
                             DE_prev_1 = 1, DE_prev_2 = 1, DE_RgivenTested_1 = 1, DE_RgivenTested_2 = 1,
                             rec_test_coverage_1 = 1, rec_test_coverage_2 = 1,
                             alpha = 0.05) {

  #checks consistency of inputs
  if(MDRI <= 0) stop("MDRI should be positive")
  if(MDRI > TIME) {warning("MDRI cannot be > TIME")
                   return(NA)}
  if(inc_2 > inc_1) {warning("inc_2 cannot be > inc_1")
                     return(NA)}
  if(any(mdrihatcov > 1, mdrihatcov < 0)) stop("check mdrihatcov range")
  if(any(frrhatcov > 1, frrhatcov < 0)) stop("check frrhatcov range")
  if(any(frrhat > 1, frrhat < 0)) stop("frrhat value must be between 0 and 1")
  if(any(inc_1 > 1, inc_1 < 0)) stop ("inc_1 value must be between 0 and 1")
  if(any(inc_2 > 1, inc_2 < 0)) stop ("inc_2 value must be between 0 and 1")
  if(any(p_pos_1 > 1, p_pos_1 < 0)) stop ("p_pos_1 value must be between 0 and 1")
  if(any(p_pos_2 > 1, p_pos_2 < 0)) stop ("p_pos_2 value must be between 0 and 1")
  if(any(rec_test_coverage_1 > 1, rec_test_coverage_1 < 0)) stop ("rec_test_coverage_1 value must be between 0 and 1")
  if(any(rec_test_coverage_2 > 1, rec_test_coverage_2 < 0)) stop ("rec_test_coverage_2 value must be between 0 and 1")

  TIME_CONV <- 365.25
  mdrihat <- MDRI/TIME_CONV
  TIME <- TIME/TIME_CONV
  p_s_1 <- 1 - p_pos_1
  p_s_2 <-1 - p_pos_2

  max_frr <- 0.1
  min_frr <- 0
  max_mdri <- 1
  min_mdri <- 0.2460
  min_n_r <- 10
  max_cov <- 0.25
  max_alpha <- 0.10
  min_prob <- 0.70

  p_r_1 <- inc_1*(1 - p_pos_1)*(mdrihat - frrhat*TIME) + frrhat*p_pos_1
  p_r_2 <- inc_2*(1 - p_pos_2)*(mdrihat - frrhat*TIME) + frrhat*p_pos_2

  p_RgivenTested_1 <- p_r_1/p_pos_1
  p_RgivenTested_2 <- p_r_2/p_pos_2

  rse_squared_1 <- data.frame()
  rse_squared_2 <- data.frame()

  rse_difference <- 1/(-qnorm(alpha/2) - qnorm(1 - power))

  rse_squared_finitesamplesize_1 <- (1/p_pos_1)*(DE_prev_1/p_s_1 + (DE_RgivenTested_1/rec_test_coverage_1)*p_RgivenTested_1*(1 - p_RgivenTested_1)/((p_RgivenTested_1 - frrhat)^2))
  rse_squared_finitesamplesize_2 <- (1/p_pos_2)*(DE_prev_2/p_s_2 + (DE_RgivenTested_2/rec_test_coverage_2)*p_RgivenTested_2*(1 - p_RgivenTested_2)/((p_RgivenTested_2 - frrhat)^2))
  rse_squared_mdri_1 <- (mdrihatcov*mdrihat/(mdrihat - frrhat*TIME))^2
  rse_squared_mdri_2 <- (mdrihatcov*mdrihat/(mdrihat - frrhat*TIME))^2
  rse_squared_frr_1 <- (frrhatcov*frrhat*(mdrihat - p_RgivenTested_1*TIME)/((mdrihat - frrhat*TIME)*(p_RgivenTested_1 - frrhat)))^2
  rse_squared_frr_2 <- (frrhatcov*frrhat*(mdrihat - p_RgivenTested_2*TIME)/((mdrihat - frrhat*TIME)*(p_RgivenTested_2 - frrhat)))^2
  lambda_1 <- ((inc_1^2)*rse_squared_finitesamplesize_1 + (inc_2^2)*rse_squared_finitesamplesize_2)

  rse_squared_1 <- c(rse_squared_finitesamplesize_1, rse_squared_mdri_1, rse_squared_frr_1)
  rse_squared_2 <- c(rse_squared_finitesamplesize_2, rse_squared_mdri_2, rse_squared_frr_2)

  cov2_1 <- sum(rse_squared_1)
  cov2_2 <- sum(rse_squared_2)

  #Estimated difference
  inc_diff <- inc_1 - inc_2

  variance_mdri <- ((mdrihatcov*mdrihat)/(mdrihat - frrhat*TIME))^2*(inc_diff^2)
  variance_frr <- ((frrhat*frrhatcov)^2/((mdrihat - frrhat*TIME)^4)*(p_pos_2*(mdrihat - p_RgivenTested_2*TIME)/p_s_2 - p_pos_1*(mdrihat - p_RgivenTested_1*TIME)/p_s_1)^2)
  min_var_diff <- variance_mdri + variance_frr
  min_rse <- sqrt(min_var_diff)/inc_diff
  max_power <- 1 - pnorm(-qnorm(alpha/2), 1/min_rse, 1)
  lambda_2 <- variance_mdri + variance_frr

  denominator <- (inc_diff^2)*(rse_difference^2) - lambda_2
  sample_size <- lambda_1 / denominator

  p_nr_1 <- 1 - p_s_1 - p_r_1
  p_nr_2 <- 1 - p_s_2 - p_r_2
  if(any(p_r_1 < 0, p_r_2 < 0, p_nr_1 < 0, p_nr_2 < 0)) {
    warning("Incompatible inputs produce negative population proportions")
    return(NA)}
  ceiling(sample_size)

}


#calculates the required power in case 2
#Case II: Assumes that the two surveys use a single MDRI estimate, but that the FRRs are independently estimated
ss_calc_case2 <- function(TIME = 730, MDRI = 200,
                          mdrihatcov = 0.05, frrhat_1 = 0.01, frrhatcov_1 = 0.20,
                          frrhat_2 = 0.01, frrhatcov_2 = 0.30,
                          inc_1 = 0.05, inc_2 = 0.03, p_pos_1 = 0.20, p_pos_2 = 0.15,
                          power = 0.8,
                          DE_prev_1 = 1, DE_prev_2 = 1, DE_RgivenTested_1 = 1, DE_RgivenTested_2 = 1,
                          rec_test_coverage_1 = 1, rec_test_coverage_2 = 1,
                          alpha = 0.05) {
  #browser()
  if(MDRI <= 0) stop("MDRI should be positive")
  if(MDRI > TIME) {warning("MDRI cannot be > TIME")
                   return(NA)}
  if(inc_2 > inc_1) {warning("inc_2 cannot be > inc_1")
                     return(NA)}
  if(any(mdrihatcov > 1, mdrihatcov < 0)) stop("check mdrihatcov range")
  if(any(frrhatcov_1 > 1, frrhatcov_1 < 0)) stop("check frrhatcov_1 range")
  if(any(frrhatcov_2 > 1, frrhatcov_2 < 0)) stop("check frrhatcov_2 range")
  if(any(frrhat_1 > 1, frrhat_1 < 0)) stop("frrhat_1 value must be between 0 and 1")
  if(any(frrhat_2 > 1, frrhat_2 < 0)) stop("frrhat_2 value must be between 0 and 1")
  if(any(inc_1 > 1, inc_1 < 0)) stop ("inc_1 value must be between 0 and 1")
  if(any(inc_2 > 1, inc_2 < 0)) stop ("inc_2 value must be between 0 and 1")
  if(any(p_pos_1 > 1, p_pos_1 < 0)) stop ("p_pos_1 value must be between 0 and 1")
  if(any(p_pos_2 > 1, p_pos_2 < 0)) stop ("p_pos_2 value must be between 0 and 1")
  if(any(rec_test_coverage_1 > 1, rec_test_coverage_1 < 0)) stop ("rec_test_coverage_1 value must be between 0 and 1")
  if(any(rec_test_coverage_2 > 1, rec_test_coverage_2 < 0)) stop ("rec_test_coverage_2 value must be between 0 and 1")

  TIME_CONV <- 365.25
  mdrihat <- MDRI/TIME_CONV
  TIME <- TIME/TIME_CONV
  p_s_1 <- 1 - p_pos_1
  p_s_2 <-1 - p_pos_2

  max_frr <- 0.1
  min_frr <- 0
  max_mdri <- 1
  min_mdri <- 0.2460
  min_n_r <- 10
  max_cov <- 0.25
  max_alpha <- 0.10
  min_prob <- 0.70

  p_r_1 <- inc_1*(1 - p_pos_1)*(mdrihat - frrhat_1*TIME) + frrhat_1*p_pos_1
  p_r_2 <- inc_2*(1 - p_pos_2)*(mdrihat - frrhat_2*TIME) + frrhat_2*p_pos_2

  p_RgivenTested_1 <- p_r_1/p_pos_1
  p_RgivenTested_2 <- p_r_2/p_pos_2

  rse_squared_1 <- data.frame()
  rse_squared_2 <- data.frame()

  rse_difference <- 1/(-qnorm(alpha/2) - qnorm(1 - power))

  rse_squared_finitesamplesize_1 <- (1/p_pos_1)*(DE_prev_1/p_s_1 + (DE_RgivenTested_1/rec_test_coverage_1)*p_RgivenTested_1*(1 - p_RgivenTested_1)/((p_RgivenTested_1 - frrhat_1)^2))
  rse_squared_finitesamplesize_2 <- (1/p_pos_2)*(DE_prev_2/p_s_2 + (DE_RgivenTested_2/rec_test_coverage_2)*p_RgivenTested_2*(1 - p_RgivenTested_2)/((p_RgivenTested_2 - frrhat_2)^2))
  rse_squared_mdri_1 <- (mdrihatcov*mdrihat/(mdrihat - frrhat_1*TIME))^2
  rse_squared_mdri_2 <- (mdrihatcov*mdrihat/(mdrihat - frrhat_2*TIME))^2
  rse_squared_frr_1 <- (frrhatcov_1*frrhat_1*(mdrihat - p_RgivenTested_1*TIME)/((mdrihat - frrhat_1*TIME)*(p_RgivenTested_1 - frrhat_1)))^2
  rse_squared_frr_2 <- (frrhatcov_2*frrhat_2*(mdrihat - p_RgivenTested_2*TIME)/((mdrihat - frrhat_2*TIME)*(p_RgivenTested_2 - frrhat_2)))^2
  lambda_1 <- ((inc_1^2)*rse_squared_finitesamplesize_1 + (inc_2^2)*rse_squared_finitesamplesize_2)

  rse_squared_1 <- c(rse_squared_finitesamplesize_1, rse_squared_mdri_1, rse_squared_frr_1)
  rse_squared_2 <- c(rse_squared_finitesamplesize_2, rse_squared_mdri_2, rse_squared_frr_2)

  cov2_1 <- sum(rse_squared_1)
  cov2_2 <- sum(rse_squared_2)

  #Estimated difference
  inc_diff <- inc_1 - inc_2

  variance_mdri<- ((mdrihat*mdrihatcov)^2)*(inc_1/(mdrihat - frrhat_1*TIME) - inc_2/(mdrihat - frrhat_2*TIME))^2
  variance_frr <- rse_squared_frr_1*inc_1^2 + rse_squared_frr_2*inc_2^2

  min_var_diff <- variance_mdri + variance_frr
  min_rse <- sqrt(min_var_diff)/inc_diff
  max_power <- 1 - pnorm(-qnorm(alpha/2), 1/min_rse, 1)
  lambda_2 <- variance_mdri + variance_frr

  denominator <- (inc_diff^2)*(rse_difference^2) - lambda_2
  sample_size <- lambda_1 / denominator

  p_nr_1 <- 1 - p_s_1 - p_r_1
  p_nr_2 <- 1 - p_s_2 - p_r_2
  if(any(p_r_1 < 0, p_r_2 < 0, p_nr_1 < 0, p_nr_2 < 0)) {
    warning("Error: Incompatible inputs produce negative population proportions")
    return(NA)}
  #sample_size
  ceiling(sample_size)

}


#calculates the required power in case 3
#Case III: Assumes that the two surveys use MDRI estimates which arise from different incidence tests, and that the FRRs are independently estimated
ss_calc_case3 <- function(TIME = 730, MDRI_1 = 200, MDRI_2 = 200,
                          mdrihatcov_1 = 0.05, mdrihatcov_2 = 0.05, frrhat_1 = 0.01, frrhatcov_1 = 0.20,
                          frrhat_2 = 0.01, frrhatcov_2 = 0.20,
                          inc_1 = 0.05, inc_2 = 0.03, p_pos_1 = 0.20, p_pos_2 = 0.15,
                          power = 0.8,
                          DE_prev_1 = 1, DE_prev_2 = 1, DE_RgivenTested_1 = 1, DE_RgivenTested_2 = 1,
                          rec_test_coverage_1 = 1, rec_test_coverage_2 = 1,
                          alpha = 0.05) {

  if(any(MDRI_1 <= 0, MDRI_2 <= 0)) stop("MDRI_1 should be positive")
  if(any(MDRI_1 > TIME, MDRI_2 > TIME)) {warning("MDRIs cannot be > TIME")
                   return(NA)}
  if(inc_2 > inc_1) {warning("inc_2 cannot be > inc_1")
                     return(NA)}
  if(any(mdrihatcov_1 > 1, mdrihatcov_1 < 0)) stop("check mdrihatcov_1 range")
  if(any(mdrihatcov_2 > 1, mdrihatcov_2 < 0)) stop("check mdrihatcov_2 range")
  if(any(frrhatcov_1 > 1, frrhatcov_1 < 0)) stop("check frrhatcov_1 range")
  if(any(frrhatcov_2 > 1, frrhatcov_2 < 0)) stop("check frrhatcov_2 range")
  if(any(frrhat_1 > 1, frrhat_1 < 0)) stop("frrhat_1 value must be between 0 and 1")
  if(any(frrhat_2 > 1, frrhat_2 < 0)) stop("frrhat_2 value must be between 0 and 1")
  if(any(inc_1 > 1, inc_1 < 0)) stop ("inc_1 value must be between 0 and 1")
  if(any(inc_2 > 1, inc_2 < 0)) stop ("inc_2 value must be between 0 and 1")
  if(any(p_pos_1 > 1, p_pos_1 < 0)) stop ("p_pos_1 value must be between 0 and 1")
  if(any(p_pos_2 > 1, p_pos_2 < 0)) stop ("p_pos_2 value must be between 0 and 1")
  if(any(rec_test_coverage_1 > 1, rec_test_coverage_1 < 0)) stop ("rec_test_coverage_1 value must be between 0 and 1")
  if(any(rec_test_coverage_2 > 1, rec_test_coverage_2 < 0)) stop ("rec_test_coverage_2 value must be between 0 and 1")

  TIME_CONV <- 365.25
  mdrihat_1 <- MDRI_1/TIME_CONV
  mdrihat_2 <- MDRI_2/TIME_CONV
  TIME <- TIME/TIME_CONV
  p_s_1 <- 1 - p_pos_1
  p_s_2 <-1 - p_pos_2

  max_frr <- 0.1
  min_frr <- 0
  max_mdri <- 1
  min_mdri <- 0.2460
  min_n_r <- 10
  max_cov <- 0.25
  max_alpha <- 0.10
  min_prob <- 0.70
#browser()
  p_r_1 <- inc_1*(1 - p_pos_1)*(mdrihat_1 - frrhat_1*TIME) + frrhat_1*p_pos_1
  p_r_2 <- inc_2*(1 - p_pos_2)*(mdrihat_2 - frrhat_2*TIME) + frrhat_2*p_pos_2

  p_RgivenTested_1 <- p_r_1/p_pos_1
  p_RgivenTested_2 <- p_r_2/p_pos_2

  rse_squared_1 <- data.frame()
  rse_squared_2 <- data.frame()

  rse_difference <- 1/(-qnorm(alpha/2) - qnorm(1 - power))

  rse_squared_finitesamplesize_1 <- (1/p_pos_1)*(DE_prev_1/p_s_1 + (DE_RgivenTested_1/rec_test_coverage_1)*p_RgivenTested_1*(1 - p_RgivenTested_1)/((p_RgivenTested_1 - frrhat_1)^2))
  rse_squared_finitesamplesize_2 <- (1/p_pos_2)*(DE_prev_2/p_s_2 + (DE_RgivenTested_2/rec_test_coverage_2)*p_RgivenTested_2*(1 - p_RgivenTested_2)/((p_RgivenTested_2 - frrhat_2)^2))
  rse_squared_mdri_1 <- (mdrihatcov_1*mdrihat_1/(mdrihat_1 - frrhat_1*TIME))^2
  rse_squared_mdri_2 <- (mdrihatcov_2*mdrihat_2/(mdrihat_2 - frrhat_2*TIME))^2
  rse_squared_frr_1 <- (frrhatcov_1*frrhat_1*(mdrihat_1 - p_RgivenTested_1*TIME)/((mdrihat_1 - frrhat_1*TIME)*(p_RgivenTested_1 - frrhat_1)))^2
  rse_squared_frr_2 <- (frrhatcov_2*frrhat_2*(mdrihat_2 - p_RgivenTested_2*TIME)/((mdrihat_2 - frrhat_2*TIME)*(p_RgivenTested_2 - frrhat_2)))^2
  lambda_1 <- ((inc_1^2)*rse_squared_finitesamplesize_1 + (inc_2^2)*rse_squared_finitesamplesize_2)

  rse_squared_1 <- c(rse_squared_finitesamplesize_1, rse_squared_mdri_1, rse_squared_frr_1)
  rse_squared_2 <- c(rse_squared_finitesamplesize_2, rse_squared_mdri_2, rse_squared_frr_2)

  cov2_1 <- sum(rse_squared_1)
  cov2_2 <- sum(rse_squared_2)

  #Estimated difference
  inc_diff <- inc_1 - inc_2
  variance_frr <- rse_squared_frr_1*inc_1^2 + rse_squared_frr_2*inc_2^2
  lambda_2 <- (inc_1^2)*(rse_squared_mdri_1 + rse_squared_frr_1) + (inc_2^2)*(rse_squared_mdri_2 + rse_squared_frr_2)

  denominator <- (inc_diff^2)*(rse_difference^2) - lambda_2
  sample_size <- lambda_1 / denominator

  p_nr_1 <- 1 - p_s_1 - p_r_1
  p_nr_2 <- 1 - p_s_2 - p_r_2
  if(any(p_r_1 < 0, p_r_2 < 0, p_nr_1 < 0, p_nr_2 < 0)) {
    warning("Error: Incompatible inputs produce negative population proportions")
    return(NA)}
  #sample_size
  ceiling(sample_size)

}


#general wrapper that takes all inputs and provides the power requested according to the case selected
#by the user, works via internal manual dispatching
#inputs and units are ad following:
#case: integer
#TIME days (corresponds to T in .xls files)
#MDRI, MDRI_1, MDRI_2 days
#mdrihatcov, mdrihatcov_1, mdrihatcov_1 covariance of MDRI estiamtes (decimal)
#frrhat, frrhat_1, frrhat_2 false recent rates, decimals <1
#frrhatcov, frrhatcov_1, frrhatcov_2 frrhat covariance (decimal)
#inc_1, inc_2, p_pos_1, p_pos_2 incidence and prevalence for the 2 survey (decimals, <1)
#n_1, n_2 number of subjects in each survey
#DE_prev_1, DE_prev_2, DE_RgivenTested_1, DE_RgivenTested_2 design effects
#rec_test_coverage_1, rec_test_coverage_2 decimal (<1)
#alpha significance level for testng (<1)
#Example of a call: power_calc(case = 3, MDRI_1 = 250, MDRI_2 = 360, frrhat_1 = 0.01, frrhat_2 = 0.04, alpha = 0.01)
#69.92957
ss_calc <- function(case = 1, TIME = 730,
                       MDRI = 200, MDRI_1 = 200, MDRI_2 = 200,
                       mdrihatcov = 0.05, mdrihatcov_1 = 0.05, mdrihatcov_2 = 0.05,
                       frrhat = 0.01, frrhat_1 = 0.01, frrhatcov = 0.20,
                      frrhatcov_1 = 0.20, frrhat_2 = 0.01, frrhatcov_2 = 0.30,
                       inc_1 = 0.05, inc_2 = 0.03, p_pos_1 = 0.20, p_pos_2 = 0.15,
                       power = 0.8,
                       DE_prev_1 = 1, DE_prev_2 = 1, DE_RgivenTested_1 = 1, DE_RgivenTested_2 = 1,
                       rec_test_coverage_1 = 1, rec_test_coverage_2 = 1,
                       alpha = 0.05) {

  #checks if case is correctly specified
  if(!sum(case == c(1, 2, 3))) {stop("Please enter a valid case value")}

  #manual dispatching according to case, passing arguments to the appropriate function
  #browser()
  if (1 == case) {
    ss <- ss_calc_case1(TIME = TIME, MDRI = MDRI,  mdrihatcov = mdrihatcov,
                        frrhat = frrhat, frrhatcov = frrhatcov,
                            inc_1 = inc_1, inc_2 = inc_2, p_pos_1 = p_pos_1, p_pos_2 = p_pos_2,
                            power = power,
                            DE_prev_1 = DE_prev_1, DE_prev_2 = DE_prev_2,
                            DE_RgivenTested_1 = DE_RgivenTested_1, DE_RgivenTested_2 = DE_RgivenTested_2,
                            rec_test_coverage_1 = rec_test_coverage_1, rec_test_coverage_2 = rec_test_coverage_2,
                            alpha = alpha)
  }

  if (2 == case) {
    ss <- ss_calc_case2(TIME = TIME, MDRI = MDRI, mdrihatcov = mdrihatcov,
                            frrhat_1 = frrhat_1, frrhatcov_1 = frrhatcov_1, frrhat_2 = frrhat_2, frrhatcov_2 = frrhatcov_2,
                            inc_1 = inc_1, inc_2 = inc_2, p_pos_1 = p_pos_1, p_pos_2 = p_pos_2,
                            power = power,
                            DE_prev_1 = DE_prev_1, DE_prev_2 = DE_prev_2,
                            DE_RgivenTested_1 = DE_RgivenTested_1, DE_RgivenTested_2 = DE_RgivenTested_2,
                            rec_test_coverage_1 = rec_test_coverage_1, rec_test_coverage_2 = rec_test_coverage_2,
                            alpha = alpha)

  }

  if(3 == case) {
    ss <- ss_calc_case3(TIME = TIME, MDRI_1 = MDRI_1, MDRI_2 = MDRI_2,
                            mdrihatcov_1 = mdrihatcov_1, mdrihatcov_2 = mdrihatcov_2,
                            frrhat_1 = frrhat_1, frrhatcov_1 = frrhatcov_1, frrhat_2 = frrhat_2, frrhatcov_2 = frrhatcov_2,
                            inc_1 = inc_1, inc_2 = inc_2, p_pos_1 = p_pos_1, p_pos_2 = p_pos_2,
                            power = power,
                            DE_prev_1 = DE_prev_1, DE_prev_2 = DE_prev_2,
                            DE_RgivenTested_1 = DE_RgivenTested_1, DE_RgivenTested_2 = DE_RgivenTested_2,
                            rec_test_coverage_1 = rec_test_coverage_1, rec_test_coverage_2 = rec_test_coverage_2,
                            alpha = alpha)

  }

  return(ss)
}
