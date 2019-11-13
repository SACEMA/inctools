# Created by and Copyright (C) 2015-2016 Stefano Ongarello (FIND).
# Recoded by Lamin Juwara (McGill) to use functions from R Package 'inctools' (2017/18)
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


# Here, we create a general wrapper that takes all the inputs and computes the sample sizes 
# for the three cases based on the biomaker test parameters (MDRI, FRR, and RSE). 
# See the BMest input in inctools.
ss_calc <- function(case = 1, TIME = 730, MDRI = 200, MDRI_1 = 200, MDRI_2 = 200,
                    mdrihatcov = 0.05, mdrihatcov_1 = 0.05, mdrihatcov_2 = 0.05, 
                    frrhat = 0.01, frrhat_1 = 0.01, frrhatcov = 0.20, frrhatcov_1 = 0.20, 
                    frrhat_2 = 0.01, frrhatcov_2 = 0.30, inc_1 = 0.05, inc_2 = 0.03,
                    p_pos_1 = 0.20, p_pos_2 = 0.15, power = 0.8, DE_prev_1 = 1, DE_prev_2 = 1,
                    DE_RgivenTested_1 = 1, DE_RgivenTested_2 = 1,
                    rec_test_coverage_1 = 1, rec_test_coverage_2 = 1, alpha = 0.05) {
  
  #checks if case is correctly specified
  if(!sum(case == c(1, 2, 3))) {stop("Please enter a valid case value")}
  
  #manual dispatching according to case, passing arguments to the appropriate function
  
  if (1 == case) { #Case I: Assumes that the two surveys use a single MDRI and FRRs estimate.
    temp_ss <- incpower( MDRI = MDRI, FRR = frrhat, BigT = TIME, RSE_FRR = frrhatcov,
                         RSE_MDRI = mdrihatcov, DE_H = c(DE_prev_1,DE_prev_2), 
                         DE_R = c(DE_RgivenTested_1,DE_RgivenTested_2), I1 = 0.05, I2 = 0.03,
                         PrevH1 =  p_pos_1, PrevH2 =  p_pos_2, Power = power, alpha = alpha,
                         CR = c(rec_test_coverage_1,rec_test_coverage_2),
                         BMest = 'same.test', SS = 'out' )
    ss<-temp_ss$Minimum.Common.SS
    
  }
  
  if (2 == case) { #Case II: Assumes that the two surveys use a single MDRI estimate, but that the FRRs are independently estimated
    temp_ss <- incpower( MDRI = MDRI, FRR = c(frrhat_1,frrhat_2), BigT = TIME, 
                         RSE_FRR = c(frrhatcov_1,frrhatcov_2),
                         RSE_MDRI = mdrihatcov, DE_H = c(DE_prev_1,DE_prev_2), 
                         DE_R = c(DE_RgivenTested_1,DE_RgivenTested_2), I1 = 0.05, I2 = 0.03,
                         PrevH1 =  p_pos_1, PrevH2 =  p_pos_2, Power = power, alpha = alpha,
                         CR = c(rec_test_coverage_1,rec_test_coverage_2),
                         BMest = 'FRR.indep', SS = 'out' )
    ss <- temp_ss$Minimum.Common.SS
      #sum(temp_ss$Implied.Subject.Counts[1:2,1])
    
  }
  
  if(3 == case) { #Case III: Assumes that the two surveys use MDRI estimates which arise from different incidence tests, and that the FRRs are independently estimated
    temp_ss <-incpower( MDRI = c(MDRI_1,MDRI_2), FRR = c(frrhat_1,frrhat_2), BigT = TIME, 
                        RSE_FRR = c(frrhatcov_1,frrhatcov_2),
                        RSE_MDRI = c(mdrihatcov_1,mdrihatcov_2), 
                        DE_H = c(DE_prev_1,DE_prev_2), 
                        DE_R = c(DE_RgivenTested_1,DE_RgivenTested_2), I1 = 0.05, I2 = 0.03,
                        PrevH1 =  p_pos_1, PrevH2 =  p_pos_2, Power = power, alpha = alpha,
                        CR = c(rec_test_coverage_1,rec_test_coverage_2),
                        BMest = 'MDRI.FRR.indep', SS = 'out' )
    ss <- temp_ss$Minimum.Common.SS
    
  }
  
  return(ss)
}

#changes the names of the input dataset to be used by ss_calc
assign_colnames <- function(cd) {
  #names are correct
  if(all(c("inc_1", "prev_1", "MDRI", "FRR", "scenario") %in% colnames(cd))) {
    colnames(cd)[which("FRR" == colnames(cd))] <- "frrhat"
    colnames(cd)[which("prev_1" == colnames(cd))] <- "p_pos_1"}
  #names are not correct but assume positions are
  else {
    colnames(cd) <- c("inc_1", "p_pos_1", "MDRI", "frrhat", "scenario")
  }
  return(cd)
}

#reshapes the df to be processed by sscalc
assemble_data <- function(ip, percent_reduction = 0.5, DE_p_1 = 1.3, DE_p_2 = 1.3, DE_R_1 = 1.3, DE_R_2 = 1.3,
                          power = 0.8, alpha = 0.05) {
  if(is.null(ip)) {stop("Please provide input")}
  #ip <- assign_colnames(ip)
  temp <- cbind(ip, inc_2 = ip$inc_1 * percent_reduction, p_pos_2 = ip$p_pos_1,
                DE_prev_1 = DE_p_1, DE_prev_2 = DE_p_2,
                DE_RgivenTested_1 = DE_R_1, DE_RgivenTested_2 = DE_R_2,
                power = power, alpha = alpha)
  return(temp)
}

#calculates the sample sizes via mdply
do_table <- function(cd, percent_reduction, DE_p_1, DE_p_2, DE_R_1, DE_R_2, power, alpha) {

  cd <- assign_colnames(cd)
  cd$inc_1 <- cd$inc_1 * 1/100
  cd$p_pos_1 <- cd$p_pos_1 * 1/100
  cd$frrhat <- cd$frrhat * 1/100
  df <- assemble_data(cd, percent_reduction = (100 - percent_reduction)/100,
                      DE_p_1 = DE_p_1, DE_p_2 = DE_p_2,
                      DE_R_1 = DE_R_1, DE_R_2 = DE_R_2,
                      power = power, alpha = alpha)
  
  res_I001 <- mdply(df[, - which(colnames(df) == "scenario")], ss_calc, case = 1)
  res <- cbind(res_I001, scenario = df[,5])
  id <- which(colnames(res) == "V1")
  colnames(res)[id] <- "N"
  #res
  res_clean <- res[, c("scenario", "inc_1", "p_pos_1", "MDRI", "frrhat", "N")]
  res_clean[, c(2, 3, 5)] <- 100*res_clean[, c(2, 3, 5)]
  res_clean$N <- round(res_clean$N, digits = 0)
  #res_clean[, 2:6] <- sapply(res_clean[, 2:6], FUN = function(x) prettyNum(x, big.mark=","))
  colnames(res_clean) <- c("Scenario", "inc_1", "prev_1", "MDRI", "FRR", "Sample size")

  res_clean_extra <- cbind(res_clean, alpha = res[, c("alpha")], power = res[, c("power")],
                           inc_2 = 100 * res[, c("inc_2")], prev_2 = 100 * res[, c("p_pos_2")],
                           #percent_reduction = res[, c("percent_reduction")],
                           DE_prev1 = res[, c("DE_prev_1")], DE_prev2 = res[, c("DE_prev_2")],
                           DE_RgivenTested1 = res[, c("DE_RgivenTested_1")], DE_RgivenTested2 = res[, c("DE_RgivenTested_2")]
                           )
  res_clean_extra$percent_reduction <- round(-100*(res_clean_extra$inc_2 - res_clean_extra$inc_1)/res_clean_extra$inc_1)
  res_clean_extra[, 2:15] <- sapply(res_clean_extra[, 2:15], FUN = function(x) prettyNum(x, big.mark=","))

  return(res_clean_extra)
}

