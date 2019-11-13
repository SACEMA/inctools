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


incdiff_calc<-function(COV_MDRI = 0, COV_FRR=0,
                       PrevH_1 = 0.20, RSE_PrevH_1 = 0.028,
                       PrevR_1 = 0.10, RSE_PrevR_1 = 0.094,
                       MDRI_1 = 200, RSE_MDRI_1 =.05,
                       FRR_1 = .01, RSE_FRR_1 = .2,
                       MDRI_2 = 180, RSE_MDRI_2 = .07,
                       FRR_2 = .009, RSE_FRR_2 = .2,
                       PrevH_2 = .21, RSE_PrevH_2 = .03,
                       PrevR_2 = .13, RSE_PrevR_2 = .095,
                       BigT = 730){

    
    #manual dispatching according to case, passing arguments to the appropriate function
    
    if (COV_MDRI==0 & COV_FRR==0) { #Case I: Assumes that the two surveys use a single MDRI and FRRs estimate.
      temp<-incprops(PrevH = c(PrevH_1,PrevH_2), 
                     RSE_PrevH = c(RSE_PrevH_1, RSE_PrevH_2),
                     PrevR = c(PrevR_1,PrevR_2), 
                     RSE_PrevR = c(RSE_PrevR_1, RSE_PrevR_2),
                     BS_Count = 10000, Boot = FALSE,
                     MDRI = c(MDRI_1,MDRI_1),
                     RSE_MDRI = c(RSE_MDRI_1,RSE_MDRI_1),
                     FRR = c(FRR_1,FRR_1),
                     RSE_FRR = c(RSE_FRR_1,RSE_FRR_1),
                     BigT = BigT, 
                     BMest = 'same.test')
      incdiff<-temp$Incidence.Difference.Statistics
    }
    
    if (COV_MDRI==0 & COV_FRR==1) { #Case II: Assumes that the two surveys use a single MDRI estimate, but FRRs are independently estimated
      temp<-incprops(PrevH = c(PrevH_1,PrevH_2), 
                     RSE_PrevH = c(RSE_PrevH_1, RSE_PrevH_2),
                     PrevR = c(PrevR_1,PrevR_2), 
                     RSE_PrevR = c(RSE_PrevR_1, RSE_PrevR_2),
                     BS_Count = 10000, Boot = FALSE,
                     MDRI = c(MDRI_1,MDRI_1),
                     RSE_MDRI = c(RSE_MDRI_1,RSE_MDRI_1),
                     FRR = c(FRR_1,FRR_2),
                     RSE_FRR = c(RSE_FRR_1,RSE_FRR_2),
                     BigT = BigT,
                     BMest = 'FRR.indep')
      incdiff <- temp$Incidence.Difference.Statistics
      
    }
  if(COV_MDRI==1 & COV_FRR==0) { #Case not allowed
    ss <- "The two surveys cannot have the same FRR and independent MDRIs"
    
  }
  
  if(COV_MDRI==1 & COV_FRR==1){ #Case III: Assumes that the two surveys use MDRI estimates which arise from different incidence tests, and that the FRRs are independently estimated
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


