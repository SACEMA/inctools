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

ss_calc <- function(COV_MDRI = 0, COV_FRR=0,  BigT = 730, 
                    MDRI_1 = 200, MDRI_2 = 200,
                    RSE_MDRI_1 = 0.05, RSE_MDRI_2 = 0.05, 
                    FRR_1 = 0.01, RSE_FRR_1 = 0.20, 
                    FRR_2 = 0.01, RSE_FRR_2 = 0.30,
                    I1 = 0.05, I2 = 0.03,
                    PrevH1 = 0.20, PrevH2 = 0.15, Power = 0.8, DE_H1 = 1, DE_H2 = 1,
                    DE_R1 = 1, DE_R2 = 1,
                    CR_1 = 1, CR_2 = 1, alpha = 0.05) {



  #manual dispatching according to case, passing arguments to the appropriate function

  if (COV_MDRI==0 & COV_FRR==0) { #Case I: Assumes that the two surveys use a single MDRI and FRRs estimate.
    temp_ss <-suppressWarnings(incpower( MDRI = MDRI_1, FRR = FRR_1, BigT = BigT, RSE_FRR = RSE_FRR_1,
                                         RSE_MDRI = RSE_MDRI_1, DE_H = c(DE_H1,DE_H2), 
                                         DE_R = c(DE_R1,DE_R2), I1 = 0.05, I2 = 0.03,
                                         PrevH1 =  PrevH1, PrevH2 =  PrevH2, Power = Power, alpha = alpha,
                                         CR = c(CR_1,CR_2),
                                         BMest = 'same.test', SS = 'out' ))  
    ss<-temp_ss$Minimum.Common.SS
  }

  if (COV_MDRI==0 & COV_FRR==1) { #Case II: Assumes that the two surveys use a single MDRI estimate, but that the FRRs are independently estimated
    temp_ss <- suppressWarnings(incpower( MDRI = MDRI_1, FRR = c(FRR_1,FRR_2), BigT = BigT, 
                                          RSE_FRR = c(RSE_FRR_1,RSE_FRR_2),
                                          RSE_MDRI = RSE_MDRI_1, DE_H = c(DE_H1,DE_H2), 
                                          DE_R = c(DE_R1,DE_R2), I1 = 0.05, I2 = 0.03,
                                          PrevH1 =  PrevH1, PrevH2 =  PrevH2, Power = Power, alpha = alpha,
                                          CR = c(CR_1,CR_2),
                                          BMest = 'FRR.indep', SS = 'out' )) 
    ss <- temp_ss$Minimum.Common.SS
    
  }
  if(COV_MDRI==1 & COV_FRR==0) { #Case not allowed
    ss <- "The two surveys cannot have the same FRR and independent MDRIs"
    
  }

  if(COV_MDRI==1 & COV_FRR==1) { #Case III: Assumes that the two surveys use MDRI estimates which arise from different incidence tests, and that the FRRs are independently estimated
    temp_ss <-suppressWarnings(incpower( MDRI = c(MDRI_1,MDRI_2), FRR = c(FRR_1,FRR_2), BigT = BigT, 
                                         RSE_FRR = c(RSE_FRR_1,RSE_FRR_2),
                                         RSE_MDRI = c(RSE_MDRI_1,RSE_MDRI_2), 
                                         DE_H = c(DE_H1,DE_H2), 
                                         DE_R = c(DE_R1,DE_R2), I1 = 0.05, I2 = 0.03,
                                         PrevH1 =  PrevH1, PrevH2 =  PrevH2, Power = Power, alpha = alpha,
                                         CR = c(CR_1,CR_2),
                                         BMest = 'MDRI.FRR.indep', SS = 'out' )) 
    ss <- temp_ss$Minimum.Common.SS

  }

  return(ss)
}

