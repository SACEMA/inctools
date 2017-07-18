# Created by and Copyright (C)  Lamin Juwara (McGill)(2017/18)
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


# We create a general wrapper that takes all the inputs and computes the sample sizes 
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
    temp_ss <-suppressWarnings(incpower( MDRI = MDRI, FRR = frrhat, BigT = TIME, RSE_FRR = frrhatcov,
                                         RSE_MDRI = mdrihatcov, DE_H = c(DE_prev_1,DE_prev_2), 
                                         DE_R = c(DE_RgivenTested_1,DE_RgivenTested_2), I1 = 0.05, I2 = 0.03,
                                         PrevH1 =  p_pos_1, PrevH2 =  p_pos_2, Power = power, alpha = alpha,
                                         CR = c(rec_test_coverage_1,rec_test_coverage_2),
                                         BMest = 'same.test', SS = 'out' ))  
    ss<-temp_ss$Minimum.Common.SS
  }

  if (2 == case) { #Case II: Assumes that the two surveys use a single MDRI estimate, but that the FRRs are independently estimated
    temp_ss <- suppressWarnings(incpower( MDRI = MDRI, FRR = c(frrhat_1,frrhat_2), BigT = TIME, 
                                          RSE_FRR = c(frrhatcov_1,frrhatcov_2),
                                          RSE_MDRI = mdrihatcov, DE_H = c(DE_prev_1,DE_prev_2), 
                                          DE_R = c(DE_RgivenTested_1,DE_RgivenTested_2), I1 = 0.05, I2 = 0.03,
                                          PrevH1 =  p_pos_1, PrevH2 =  p_pos_2, Power = power, alpha = alpha,
                                          CR = c(rec_test_coverage_1,rec_test_coverage_2),
                                          BMest = 'FRR.indep', SS = 'out' )) 
    ss <- temp_ss$Minimum.Common.SS
    
  }

  if(3 == case) { #Case III: Assumes that the two surveys use MDRI estimates which arise from different incidence tests, and that the FRRs are independently estimated
    temp_ss <-suppressWarnings(incpower( MDRI = c(MDRI_1,MDRI_2), FRR = c(frrhat_1,frrhat_2), BigT = TIME, 
                                         RSE_FRR = c(frrhatcov_1,frrhatcov_2),
                                         RSE_MDRI = c(mdrihatcov_1,mdrihatcov_2), 
                                         DE_H = c(DE_prev_1,DE_prev_2), 
                                         DE_R = c(DE_RgivenTested_1,DE_RgivenTested_2), I1 = 0.05, I2 = 0.03,
                                         PrevH1 =  p_pos_1, PrevH2 =  p_pos_2, Power = power, alpha = alpha,
                                         CR = c(rec_test_coverage_1,rec_test_coverage_2),
                                         BMest = 'MDRI.FRR.indep', SS = 'out' )) 
    ss <- temp_ss$Minimum.Common.SS

  }

  return(ss)
}
