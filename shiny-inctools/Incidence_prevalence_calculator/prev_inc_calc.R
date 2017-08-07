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

# prevalence calculation function based on the inctools function prevcount
prevalence_calc <- function(N = 5000, N_H = 1000,
                            N_testR = 1000, N_R = 50,
                            DE_H = 1, DE_R = 1,
                            # Boot = TRUE,
                            # BMest = 'same.test',
                            MDRI = 210, RSE_MDRI = 0.05, FRR = 0.005, RSE_FRR = 0.19,
                            BigT = 700){
  
  temp<-prevcounts(N = N, N_H = N_H, N_testR = N_testR, N_R = N_R, 
                   DE_H = DE_H, DE_R = DE_R)
  return(temp)
}

# Instantaneous incidence calculator based on the inctools function inccounts
incidence_calc <- function(N = 5000, N_H = 1000,
                            N_testR = 1000, N_R = 50,
                            DE_H = 1, DE_R = 1,
                            # Boot = FALSE,
                            # BMest = 'same.test',
                            MDRI = 210, RSE_MDRI = 0.05, FRR = 0.005, RSE_FRR = 0.19,
                            BigT = 700){
  
  temp<-inccounts(N = c(N), N_H = N_H,
                  N_testR = N_testR, N_R = N_R,
                  DE_H = DE_H, DE_R = DE_R,
                  Boot = FALSE,
                  BMest = 'same.test',
                  MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR,
                  BigT = BigT)
  return(temp$Incidence.Statistics)
}

# Annual Risk of infection calculator based on inccounts
risk_of_infection_calc <- function(N = 5000, N_H = 1000,
                           N_testR = 1000, N_R = 50,
                           DE_H = 1, DE_R = 1,
                           # Boot = TRUE,
                           # BMest = 'same.test',
                           MDRI = 210, RSE_MDRI = 0.05, 
                           FRR = 0.005, RSE_FRR = 0.19,
                           BigT = 700){
  
  temp<-inccounts(N = c(N), N_H = N_H,
                  N_testR = N_testR, N_R = N_R,
                  DE_H = DE_H, DE_R = DE_R,
                  Boot = TRUE,
                  BMest = 'same.test',
                  MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR,
                  BigT = BigT)
  return(temp$Annual.Risk.of.Infection)
}

# A general wrapper/function making use of the incidence, prevalence and risk functions from above
prev_inc_calc_counts<-function(N = 5000, N_H = 1000,
                               N_testR = 1000, N_R = 50,
                               DE_H = 1, DE_R = 1,
                               MDRI = 210, RSE_MDRI = 0.05,
                               FRR = 0.005, RSE_FRR = 0.19,
                               BigT = 700){
  temp<-cbind(
  prevalence_calc(N = N, N_H = N_H, N_testR = N_testR , N_R = N_R,
                  DE_H = DE_H, DE_R = DE_R, MDRI = MDRI, RSE_MDRI = RSE_MDRI,
                  FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT),
  incidence_calc(N = N, N_H = N_H, N_testR = N_testR , N_R = N_R,
                 DE_H = DE_H, DE_R = DE_R, MDRI = MDRI, RSE_MDRI = RSE_MDRI,
                 FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT),
  risk_of_infection_calc(N = N, N_H = N_H, N_testR = N_testR , N_R = N_R,
                         DE_H = DE_H, DE_R = DE_R, MDRI = MDRI, RSE_MDRI = RSE_MDRI,
                         FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT)
  )
  return(temp)
}


# calculator for incidence and annual risk of infection using the
# function inprops and proportions as inputs.
prev_inc_calc_incprop <- function(PrevH = 0.20, RSE_PrevH = 0.028,
                                  PrevR = 0.10, RSE_PrevR = 0.09,
                                  MDRI = 210, RSE_MDRI = 0.05,
                                  FRR = 0.005, RSE_FRR = 0.19,
                                  BigT = 700){
  temp<-incprops(PrevH = PrevH, RSE_PrevH = RSE_PrevH,
                 PrevR = PrevR, RSE_PrevR = RSE_PrevR,
                 BS_Count = 10000, Boot = TRUE, MDRI = MDRI,
                 RSE_MDRI = RSE_MDRI, FRR = FRR,
                 RSE_FRR = RSE_FRR, BigT = BigT)
  df.values<-cbind(temp$Incidence.Statistics,temp$Annual.Risk.of.Infection)
  return(df.values)
}


