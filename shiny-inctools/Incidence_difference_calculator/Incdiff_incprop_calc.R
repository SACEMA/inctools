# Created by and Copyright (C) Lamin Juwara (McGill)(2017/18)
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

# Incidence difference estimate function


incdiff_calc<-function(case = 1, survey_number = 2,
                       PrevH_1 = 0.20, RSE_PrevH_1 = 0.028,
                       PrevR_1 = 0.10, RSE_PrevR_1 = 0.094,
                       MDRI = 200, RSE_MDRI =.05,
                       FRR = .01, RSE_FRR = .2,
                       MDRI_1 = 200, RSE_MDRI_1 =.05,
                       FRR_1 = .01, RSE_FRR_1 = .2,
                       PrevH_2 = .21, RSE_PrevH_2 = .03,
                       PrevR_2 = .13, RSE_PrevR_2 = .095,
                       MDRI_2 = 180, RSE_MDRI_2 = .07,
                       FRR_2 = .009, RSE_FRR_2 = .2,
                       BigT = 730){
    # checks if case is correctly specified
    if(!sum(case == c(1, 2, 3))) {stop("Please enter a valid case value")}
    
    #manual dispatching according to case, passing arguments to the appropriate function
    
    if (1 == case) { #Case I: Assumes that the two surveys use a single MDRI and FRRs estimate.
      temp<-incprops(PrevH = c(PrevH_1,PrevH_2), 
                     RSE_PrevH = c(RSE_PrevH_1, RSE_PrevH_2),
                     PrevR = c(PrevR_1,PrevR_2), 
                     RSE_PrevR = c(RSE_PrevR_1, RSE_PrevR_2),
                     BS_Count = 10000, Boot = FALSE,
                     MDRI = c(MDRI,MDRI),
                     RSE_MDRI = c(RSE_MDRI,RSE_MDRI),
                     FRR = c(FRR,FRR),
                     RSE_FRR = c(RSE_FRR,RSE_FRR),
                     BigT = BigT, 
                     BMest = 'same.test')
      incdiff<-temp$Incidence.Difference.Statistics
    }
    
    if (2 == case) { #Case II: Assumes that the two surveys use a single MDRI estimate, but that the FRRs are independently estimated
      temp<-incprops(PrevH = c(PrevH_1,PrevH_2), 
                     RSE_PrevH = c(RSE_PrevH_1, RSE_PrevH_2),
                     PrevR = c(PrevR_1,PrevR_2), 
                     RSE_PrevR = c(RSE_PrevR_1, RSE_PrevR_2),
                     BS_Count = 10000, Boot = FALSE,
                     MDRI = c(MDRI,MDRI),
                     RSE_MDRI = c(RSE_MDRI,RSE_MDRI),
                     FRR = c(FRR_1,FRR_2),
                     RSE_FRR = c(RSE_FRR_1,RSE_FRR_2),
                     BigT = BigT,
                     BMest = 'FRR.indep')
      incdiff <- temp$Incidence.Difference.Statistics
      
    }
    
    if(3 == case) { #Case III: Assumes that the two surveys use MDRI estimates which arise from different incidence tests, and that the FRRs are independently estimated
      temp<-incprops(PrevH = c(PrevH_1,PrevH_2), 
                     RSE_PrevH = c(RSE_PrevH_1, RSE_PrevH_2),
                     PrevR = c(PrevR_1,PrevR_2), 
                     RSE_PrevR = c(RSE_PrevR_1, RSE_PrevR_2),
                     BS_Count = 10000, Boot = FALSE,
                     MDRI = c(MDRI_1,MDRI_2),
                     RSE_MDRI = c(RSE_MDRI_1,RSE_MDRI_2),
                     FRR = c(FRR_1,FRR_2),
                     RSE_FRR = c(RSE_FRR_1,RSE_FRR_2),
                     BigT = BigT,
                     BMest = 'MDRI.FRR.indep')
      incdiff <- temp$Incidence.Difference.Statistics
      
    }
  

    return(incdiff)
}


