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

#calculates the incidence difference and its cov from two surveys based on pupolation size and assay characteristics
inc_diff <- function(MDRI = 200, TIME = 730, frrhat = 0.01,
                     mdrihatcov = 0.05, frrhatcov = 0.2,
                     n_s_1 = 4000, n_s_2 = 4100,
                     n_t_1 = 1000, n_t_2 = 950,
                     n_pos_1 = 1000, n_pos_2 = 950,
                     n_r_1 = 70, n_r_2 = 45,
                     DE_prev_1 = 1, DE_RgivenTested_1 = 1,
                     DE_prev_2 = 1, DE_RgivenTested_2 = 1) {
  mdrihat <- MDRI/365.25
  TIME <- TIME /365.25
  n_1 <- n_s_1 + n_pos_1
  n_2 <- n_s_2 + n_pos_2
  #HIV-negative (proportion)
  p_s_1 <- n_s_1/n_1
  p_s_2 <- n_s_2/n_2
  
  # temp_inc_diff<- incprops(PrevH = c(0.20,0.21), RSE_PrevH = c(0.028,0.03),
  #                          PrevR = c(0.10,0.13), RSE_PrevR = c(0.094,0.095),
  #                          BS_Count = 10000, Boot = FALSE,
  #                          MDRI = 200, RSE_MDRI = 0.05, 
  #                          FRR = 0.01, RSE_FRR = 0.2, BigT = 730)
  # temp_inc_diff$Incidence.Difference.Statistics[1,c("Diff","RSE.Diff")]
  

  mdrihat <- MDRI/365.25
  TIME <- TIME /365.25
  n_1 <- n_s_1 + n_pos_1
  n_2 <- n_s_2 + n_pos_2
  #HIV-negative (proportion)
  p_s_1 <- n_s_1/n_1
  p_s_2 <- n_s_2/n_2
  #HIV-positive and "recent"
  p_r_1 <- n_r_1/n_1
  p_r_2 <- n_r_2/n_2
  #HIV-positive and "not recent"
  j52 <- (n_t_1 - n_r_1)/n_1
  k52 <- (n_t_2 - n_r_2)/n_2
  #HIV-positive and not tested for recency
  j53 <- (n_pos_1 - n_t_1)/n_1
  k53 <- (n_pos_2 - n_t_2)/n_2
  #HIV positive
  p_pos_1 <- n_pos_1/n_1
  p_pos_2 <- n_pos_2/n_2
  p_pool <- (p_pos_1*n_1 + p_pos_2*n_2)/(n_1 + n_2)
  #HIV-'recent'|(HIV-postive and tested for recency)
  p_RgivenTested_1 <- n_r_1/n_t_1
  p_RgivenTested_2 <- n_r_2/n_t_2

  #prevalence standard error
  prev_se_1 <- sqrt(p_pos_1*(1 - p_pos_1)/n_1)
  prev_se_2 <- sqrt(p_pos_2*(1 - p_pos_2)/n_2)

  inc_num_1 <- p_pos_1*(p_RgivenTested_1 - frrhat)
  inc_num_2 <- p_pos_2*(p_RgivenTested_2 - frrhat)
  inc_den_1 <- p_s_1*(mdrihat - TIME*frrhat)
  inc_den_2 <- p_s_2*(mdrihat - TIME*frrhat)
  j63 <- inc_num_1/inc_den_1
  k63 <- inc_num_2/inc_den_2

  #Components of incidence estimator RSE
  j65 <- DE_prev_1/(n_1*p_pos_1*(1 - p_pos_1)) + (DE_RgivenTested_1*p_RgivenTested_1*(1 - p_RgivenTested_1))/(n_t_1*(p_RgivenTested_1 - frrhat)^2)
  k65 <- DE_prev_2/(n_2*p_pos_2*(1 - p_pos_2)) + (DE_RgivenTested_2*p_RgivenTested_2*(1 - p_RgivenTested_2))/(n_t_2*(p_RgivenTested_2 - frrhat)^2)
  j66 <- (mdrihatcov*mdrihat/(mdrihat - frrhat*TIME))^2
  k66 <- j66
  j67 <- (frrhatcov*frrhat*(mdrihat - p_RgivenTested_1*TIME)/((mdrihat - frrhat*TIME)*(p_RgivenTested_1 - frrhat)))^2
  k67 <- (frrhatcov*frrhat*(mdrihat - p_RgivenTested_2*TIME)/((mdrihat - frrhat*TIME)*(p_RgivenTested_2 - frrhat)))^2
  cov2_1 <- j65 + j66 + j67
  cov2_2 <- k65 + k66 + k67
  u_inc1 <- p_pos_1*(p_RgivenTested_1 - frrhat)/(p_s_1*(mdrihat - frrhat*TIME))
  u_inc2 <- p_pos_2*(p_RgivenTested_2 - frrhat)/(p_s_2*(mdrihat - frrhat*TIME))
  #Components of difference estimator RSE
  n65 <- (u_inc1^2*j65 + u_inc2^2*k65)
  n66 <- (((mdrihatcov*mdrihat)/(mdrihat - frrhat*TIME))^2)*(u_inc1 - u_inc2)^2
  n67 <- (((frrhatcov*frrhat)^2)/(mdrihat - frrhat*TIME)^4)*(((p_pos_1*(mdrihat - p_RgivenTested_1*TIME)/p_s_1)-(p_pos_2*(mdrihat - p_RgivenTested_2*TIME)/p_s_2))^2)
  #Incidence difference
  inc_diff <- u_inc1 - u_inc2
  #Incidence difference covariance
  inc_diff_cov <- sqrt((n65 + n66 + n67))/abs(inc_diff)
  return(data.frame(inc_diff, inc_diff_cov))
}


inc_diff_inverted <- function(MDRI = 200, TIME = 730, frrhat = 0.01,
                              mdrihatcov = 0.05, frrhatcov = 0.2,
                              inc_1 = 0.015, inc_2 = 0.012, 
                              prev_1 = 0.15, prev_2 = 0.11, inc_diff_cov = 0.3,
                              cov2_1 = 0.05, cov2_2 = 0.05,
                              DE_prev_1 = 1, DE_RgivenTested_1 = 1, 
                              DE_prev_2 = 1, DE_RgivenTested_2 = 1) {

  inc_diff <- inc_2 - inc_1
  mdrihat <- MDRI/365.25
  TIME <- TIME /365.25
  #proportion HIV negative
  p_neg_1 <- 1 - prev_1
  p_neg_2 <- 1 - prev_2
  #proportion HIV positive and recent
  p_r_1 <- inc_1 * (1 - prev_1) * (mdrihat - frrhat*TIME) + prev_1*frrhat
  p_r_2 <- inc_2 * (1 - prev_2) * (mdrihat - frrhat*TIME) + prev_2*frrhat
  #proportion of HIV positive
  p_pos_1 <- prev_1
  p_pos_2 <- prev_2
#  p_pos_1 <- 1 - p_neg_1 - p_r_1
#  p_pos_2 <- 1 - p_neg_2 - p_r_2
  #proportions of recently infected given being tested
  p_RgivenTested_1 <- p_r_1/p_pos_1
  p_RgivenTested_2 <- p_r_2/p_pos_2
  #proportion of HIV positives and not recent
  j52 <- 1 - p_neg_1 - p_r_1
  k52 <- 1 - p_neg_2 - p_r_2

  p_t_1 <- p_pos_1 * 1
  p_t_2 <- p_pos_2 * 1


  j66 <- (mdrihatcov*mdrihat/(mdrihat - frrhat*TIME))^2
  k66 <- j66
  j67 <- (frrhatcov*frrhat*(mdrihat - p_RgivenTested_1*TIME)/((mdrihat - frrhat*TIME)*(p_RgivenTested_1 - frrhat)))^2
  k67 <- (frrhatcov*frrhat*(mdrihat - p_RgivenTested_2*TIME)/((mdrihat - frrhat*TIME)*(p_RgivenTested_2 - frrhat)))^2
  j65 <- cov2_1 -j66 - j67
  k65 <- cov2_2 - k66 - k67


  j65 <- DE_prev_1/(n_1*p_pos_1*(1 - p_pos_1)) + (DE_RgivenTested_1*p_RgivenTested_1*(1 - p_RgivenTested_1))/(n_t_1*(p_RgivenTested_1 - frrhat)^2)

  n_1 <- (1/j65) / ((DE_prev_1/(p_pos_1*(1 - p_pos_1)))  + ((DE_RgivenTested_1*p_RgivenTested_1*(1 - p_RgivenTested_1))/ ((p_RgivenTested_1 - frrhat)^2)))
  n_2 <- (1/k65) / ((DE_prev_2/(p_pos_2*(1 - p_pos_2)))  + ((DE_RgivenTested_2*p_RgivenTested_2*(1 - p_RgivenTested_2))/ ((p_RgivenTested_2 - frrhat)^2)))


  j65 <- DE_prev_1/(n_1*p_pos_1*(1 - p_pos_1)) + (DE_RgivenTested_1*p_RgivenTested_1*(1 - p_RgivenTested_1))/(p_t_1*n_1*(p_RgivenTested_1 - frrhat)^2)

  k65 <- DE_prev_2/(n_2*p_pos_2*(1 - p_pos_2)) + (DE_RgivenTested_2*p_RgivenTested_2*(1 - p_RgivenTested_2))/(p_t_2*n_2*(p_RgivenTested_2 - frrhat)^2)

}

#Calculates the incidence estimator relative standard error (RSE) implied by test characteristics (and their uncertainty), in a chosen context
test_performance <- function(MDRI = 200, TIME = 730, frrhat = 0.01,
                             mdrihatcov = 0.05, frrhatcov = 0.2,
                             I = 0.015, P = 0.15, n = 5000,
                             DE_prev = 1, DE_R = 1, rec_cov = 1,
                             alpha = 0.05) {

  TIME <- TIME/365.25
  mdrihat <- MDRI/365.25
  z <- -qnorm(alpha/2)
  p_s <- 1 - P
  p_r <- I*(1-P)*(mdrihat - frrhat*TIME) + frrhat*P
  p_nr <- 1 - p_r - p_s
  p_pos <- p_r + p_nr
  p_RgivenTested <- p_r/p_pos

  n_t <- n*p_pos*rec_cov
  cov2_n <- DE_prev/(n*p_pos*(1 - p_pos)) + (DE_R*p_RgivenTested*(1 - p_RgivenTested))/(n_t*(p_RgivenTested - frrhat)^2)
  cov2_mdri <- (mdrihatcov*mdrihat/(mdrihat - frrhat*TIME))^2
  cov2_frr <- (frrhatcov*frrhat*(mdrihat - p_RgivenTested*TIME)/((mdrihat - frrhat*TIME)*(p_RgivenTested - frrhat)))^2
  u_cov2 <- cov2_n + cov2_mdri + cov2_frr

  inc_cov <- sqrt(u_cov2)
  low_CI <- I - z*I*inc_cov
  upp_CI <- I + z*inc_cov*I
  return(data.frame(inc_cov, low_CI, upp_CI))
}

#inverts the previous function and returns sample size based on desired cov
test_performance_inverted <- function(MDRI = 180, TIME = 730, frrhat = 0.01,
                             mdrihatcov = 0.05, frrhatcov = 0.5,
                             inc = 0.015, p_pos = 0.15, inc_cov = 0.25,
                             DE_prev = 1, DE_RgivenTested = 1, rec_cov = 1,
                             alpha = 0.05) {

  TIME <- TIME/365.25
  mdrihat <- MDRI/365.25
  z <- -qnorm(alpha/2)
  p_s <- 1 - p_pos
  p_r <- inc*(1-p_pos)*(mdrihat - frrhat*TIME) + frrhat*p_pos
  p_nr <- 1 - p_r - p_s
  p_pos <- p_r + p_nr
  p_RgivenTested <- p_r/p_pos

  cov2_mdri <- (mdrihatcov*mdrihat/(mdrihat - frrhat*TIME))^2
  cov2_frr <- (frrhatcov*frrhat*(mdrihat - p_RgivenTested*TIME)/((mdrihat - frrhat*TIME)*(p_RgivenTested - frrhat)))^2
  u_cov2 <- inc_cov^2
  cov2_n <- u_cov2 - cov2_mdri - cov2_frr

  n <- (1/cov2_n) * (DE_prev/(p_pos*(1 - p_pos)) + (DE_RgivenTested*p_RgivenTested*(1 - p_RgivenTested))/((p_pos*rec_cov)*(p_RgivenTested - frrhat)^2))
  n_t <- n*p_pos*rec_cov
  low_CI <- inc - z*inc*inc_cov
  upp_CI <- inc + z*inc_cov*inc
  if(n  < 0) n <- NA
  return(data.frame(n = ceiling(n), inc_cov, low_CI, upp_CI))
}

#changes the names of the input dataset to be used by ss_calc
assign_colnames <- function(cd) {
  #names are correct
  if(all(c("inc", "prev", "MDRI", "FRR", "scenario") %in% colnames(cd))) {
    colnames(cd)[which("FRR" == colnames(cd))] <- "frrhat"
    colnames(cd)[which("prev" == colnames(cd))] <- "p_pos"}
  #names are not correct but assume positions are
  else {
    colnames(cd) <- c("inc", "p_pos", "MDRI", "frrhat", "scenario")
  }
  return(cd)
}

#reshapes the df to be processed by sscalc
assemble_data <- function(ip, DE_prev = 1.3, DE_RgivenTested = 1.3,
                          alpha = 0.05, inc_cov = 0.25, rec_cov = 1, TIME = 730,
                          mdrihatcov = 0.05, frrhatcov = 0.5) {
  if(is.null(ip)) {stop("Please provide input")}
  #ip <- assign_colnames(ip)
  temp <- cbind(ip,
                DE_prev = DE_prev, DE_RgivenTested = DE_RgivenTested,
                 alpha = alpha, inc_cov = inc_cov, rec_cov = rec_cov, TIME = TIME,
                mdrihatcov = mdrihatcov, frrhatcov = frrhatcov)
  return(temp)
}

#calculates the sample sizes via mdply
do_table <- function(cd, DE_prev, DE_RgivenTested, alpha, inc_cov, rec_cov, TIME, mdrihatcov, frrhatcov) {
  #browser()
  cd <- assign_colnames(cd)
  cd$inc <- cd$inc * 1/100
  cd$p_pos <- cd$p_pos * 1/100
  cd$frrhat <- cd$frrhat * 1/100

  df <- assemble_data(cd, DE_prev = DE_prev, DE_RgivenTested = DE_RgivenTested,
                      alpha = alpha, inc_cov = inc_cov, rec_cov = rec_cov, TIME = TIME,
                      mdrihatcov = mdrihatcov, frrhatcov = frrhatcov)

  res_I001 <- mdply(df[, - which(colnames(df) == "scenario")], test_performance_inverted)

  res <- cbind(res_I001, scenario = df[, 3])
  id <- which(colnames(res) == "n")
  colnames(res)[id] <- "N"
  id <- which(colnames(res) == "inc_cov")
  colnames(res)[id] <- "inc_RSE"
  #res
  res_clean <- res[, c("scenario", "inc", "p_pos", "MDRI", "frrhat", "N", "inc_RSE", "low_CI", "upp_CI")]
  res_clean[, c("inc", "p_pos", "frrhat", "inc_RSE", "low_CI", "upp_CI")] <- 100*res_clean[, c("inc", "p_pos", "frrhat", "inc_RSE", "low_CI", "upp_CI")]
  res_clean$N <- round(res_clean$N, digits = 0)

  res_clean_extra <- cbind(res_clean, alpha = res[, c("alpha")],
                           DE_prev = res[, c("DE_prev")], DE_RgivenTested = res[, c("DE_RgivenTested")])


  res_clean_extra[, 2:12] <- sapply(res_clean_extra[, 2:12], FUN = function(x) prettyNum(x, big.mark=","))
  colnames(res_clean_extra)[c(3,5)] <- c("prev", "FRR")

  return(res_clean_extra)
}
