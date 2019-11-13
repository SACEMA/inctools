# Created by and Copyright (C) 2015-2016 Stefano Ongarello (FIND). 
# Recoded by Lamin Juwara (McGill)(2017/18)
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

ss_calc <- function(case = 1, BigT = 730, MDRI = 200, MDRI_1 = 200, MDRI_2 = 200,
                    RSE_MDRI = 0.05, RSE_MDRI_1 = 0.05, RSE_MDRI_2 = 0.05, 
                    FRR = 0.01, FRR_1 = 0.01, RSE_FRR = 0.20, RSE_FRR_1 = 0.20, 
                    FRR_2 = 0.01, RSE_FRR_2 = 0.30, I_1 = 0.05, I_2 = 0.03,
                    PrevH_1 = 0.20, PrevH_2 = 0.15, Power = 0.8, DE_H_1 = 1, DE_H_2 = 1,
                    DE_R_1 = 1, DE_R_2 = 1,
                    CR_1 = 1, CR_2 = 1, alpha = 0.05) {

  #checks if case is correctly specified
  if(!sum(case == c(1, 2, 3))) {stop("Please enter a valid case value")}

  #manual dispatching according to case, passing arguments to the appropriate function

  if (1 == case) { #Case I: Assumes that the two surveys use a single MDRI and FRRs estimate.
    temp_ss <- incpower( MDRI = MDRI, FRR = FRR, BigT = BigT, RSE_FRR = RSE_FRR,
                    RSE_MDRI = RSE_MDRI, DE_H = c(DE_H_1,DE_H_2), 
                    DE_R = c(DE_R_1,DE_R_2), I1 = 0.05, I2 = 0.03,
                    PrevH1 =  PrevH_1, PrevH2 =  PrevH_2, Power = Power, alpha = alpha,
                    CR = c(CR_1,CR_2),
                    BMest = 'same.test', SS = 'out' )
    ss<-temp_ss$Minimum.Common.SS
#  ss=sum(temp_ss$Implied.Subject.Counts[1:2,1])
  }

  if (2 == case) { #Case II: Assumes that the two surveys use a single MDRI estimate, but that the FRRs are independently estimated
    temp_ss <- incpower( MDRI = MDRI, FRR = c(FRR_1,FRR_2), BigT = BigT, 
                         RSE_FRR = c(RSE_FRR_1,RSE_FRR_2),
                         RSE_MDRI = RSE_MDRI, DE_H = c(DE_H_1,DE_H_2), 
                         DE_R = c(DE_R_1,DE_R_2), I1 = 0.05, I2 = 0.03,
                         PrevH1 =  PrevH_1, PrevH2 =  PrevH_2, Power = Power, alpha = alpha,
                         CR = c(CR_1,CR_2),
                         BMest = 'FRR.indep', SS = 'out' )
    ss <- temp_ss$Minimum.Common.SS
    
  }

  if(3 == case) { #Case III: Assumes that the two surveys use MDRI estimates which arise from different incidence tests, and that the FRRs are independently estimated
    temp_ss <-incpower( MDRI = c(MDRI_1,MDRI_2), FRR = c(FRR_1,FRR_2), BigT = BigT, 
                                  RSE_FRR = c(RSE_FRR_1,RSE_FRR_2),
                                  RSE_MDRI = c(RSE_MDRI_1,RSE_MDRI_2), 
                                  DE_H = c(DE_H_1,DE_H_2), 
                                  DE_R = c(DE_R_1,DE_R_2), I1 = 0.05, I2 = 0.03,
                                  PrevH1 =  PrevH_1, PrevH2 =  PrevH_2, Power = Power, alpha = alpha,
                                  CR = c(CR_1,CR_2),
                                  BMest = 'MDRI.FRR.indep', SS = 'out' )
    ss <- temp_ss$Minimum.Common.SS

  }

  return(ss)
}
