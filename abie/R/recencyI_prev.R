# Assay-Based Incidence Estimation
# Copyright (C) 2015-2016, DST/NRF Centre of Excellence in Epidemiological Modelling and Analysis (SACEMA)
# and individual contributors.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

####== Incidence and incidence difference from Trinomial Counts / Prevalences of HIV and recency
#######################################################################################################
####== N:        sample size of survey 1&2 (vector)
####== N_H:      number of HIV positive of survey 1&2 (vector)
####== N_testR:  number tested for recency of survey 1&2 (vector)
####== N_R:      number of recent cases  of survey 1&2 (vector)
####== Covar_HR: covariance of probs being postive & being recent of survey 1&2 (vector)
####== MDRI:     mean duration of recent infection [days] (vector/integer)
####== RSE_MDRI: Relative standard error of MDRI   [days] (vector/integer)
####== FRR:      False recent rate (vector/integer)
####== RSE_FRR:  Relative standard error of FRR (vector/integer)
####== BigT:     post-infection time cut-off true vs false recent [days] default 730 days (integer)
####== DE_H:     Design effect of HIV-prevalence test (vector/integer)
####== DE_R:     Design effect of recency test (vector/integer)
####== BS_Vars:  Variables to be bootstrapped e.g.c("PrevH", "PrevR", "MDRI", "FRR") (string vector)
####== BS_Count: Number of iterations for bootstrapping (integer)
####== BMest:    Biomarker estimation by one the 3 options "same.test"(=default), "FRR.indep", "MDRI.FRR.idep" (string)
#######################################################################################################



prevBYcounts <- function (DE_H=1, DE_R=1, N, N_H, N_testR, N_R) {


  #######################################################################################################
  no_s <-length(N)
  if (length(DE_H)==1)     {DE_H <- rep(DE_H, times=no_s)}         else {DE_H=DE_H}
  if (length(DE_R)==1)     {DE_R <- rep(DE_R, times=no_s)}         else {DE_R=DE_R}
  stopifnot(no_s==length(N_H)  & no_s==length(N_R) & no_s==length(N_testR) &
            no_s==length(DE_H) & no_s==length(DE_R))
  PrevH <- N_H/N
  PrevR <- N_R/N_testR
  RSE_PrevH <- sqrt(((PrevH*(1-PrevH))/N)*DE_H)/PrevH
  RSE_PrevR <- sqrt(((PrevR*(1-PrevR))/N_testR)*DE_R)/PrevR

 output <- data.frame("PrevH"=PrevH, "PrevR"=PrevR, "RSE_PrevH"=RSE_PrevH, "RSE_PrevR"=RSE_PrevR)
 return(output)
}
#######################################################################################################

I_EST <- function (prevH, prevR, mdri, frr, bigt){
  prevH*(prevR - frr)/((1 - prevH)*(mdri-frr*bigt))
}

BS_SURVEY_ESTS <- function (prevH, prevR, mdri, frr, bs_count,
                            bs_var_prevH, bs_var_prevR, bs_var_mdri, bs_var_frr, covar_HR) {
  Mu    <- c(prevH, prevR, mdri, frr)
  sigma <- matrix( c(bs_var_prevH,  covar_HR,     0,           0,
                     covar_HR,      bs_var_prevR, 0,           0,
                     0,             0,            bs_var_mdri, 0,
                     0,             0,            0,           bs_var_frr),
                 nrow=4, ncol=4)
  BS_RootEst <- MASS::mvrnorm(n=bs_count, mu=Mu, Sigma=sigma, empirical=TRUE)
  return (BS_RootEst)
}

DM_FirstOrderTerms <-function (prevH, prevR, mdri, frr, bigt)   {
  fot_prevH <- (prevR-frr)/(((1-prevH)^2)*(mdri-frr*bigt))
  fot_prevR <- prevH/((1-prevH)*(mdri-frr*bigt))
  fot_mdri  <- (frr*prevH-prevR*prevH)/((1-prevH)*((mdri-frr*bigt)^2))
  fot_frr   <- (prevH*(bigt*prevR-mdri))/((1-prevH)*((mdri-frr*bigt)^2))
  return (c(fot_prevH, fot_prevR, fot_mdri, fot_frr))
}

DM_VAR_deltaI <- function (bmest, fot_prevH1, fot_prevR1, fot_mdri1, fot_frr1,
                           fot_prevH2, fot_prevR2, fot_mdri2, fot_frr2,
                           dm_var_prevH1, dm_var_prevR1, dm_var_mdri1, dm_var_frr1,
                           dm_var_prevH2, dm_var_prevR2, dm_var_mdri2, dm_var_frr2) {
  if(bmest=="sameTest"){
    DM_Var_deltaI <- ((fot_prevH1^2)*dm_var_prevH1) + ((fot_prevH2^2)*dm_var_prevH2) +
      ((fot_prevR1^2)*dm_var_prevR1) + ((fot_prevR2^2)*dm_var_prevR2) +
      ((fot_mdri1-fot_mdri2)^2*dm_var_mdri1) +
      ((fot_frr1-fot_frr2)^2*dm_var_frr1)
  }
  else if(bmest=="FRR.indep") {
    DM_Var_deltaI <- ((fot_prevH1^2)*dm_var_prevH1) + ((fot_prevH2^2)*dm_var_prevH2) +
      ((fot_prevR1^2)*dm_var_prevR1) + ((fot_prevR2^2)*dm_var_prevR2) +
      ((fot_mdri1-fot_mdri2)^2*dm_var_mdri1) +
      ((fot_frr1^2)*dm_var_frr1)     + ((fot_frr2^2)*dm_var_frr2)
  }
  else if(bmest=="MDRI.FRR.indep") {
    DM_Var_deltaI <- ((fot_prevH1^2)*dm_var_prevH1) + ((fot_prevH2^2)*dm_var_prevH2) +
      ((fot_prevR1^2)*dm_var_prevR1) + ((fot_prevR2^2)*dm_var_prevR2) +
      ((fot_mdri1^2)*dm_var_mdri1)   + ((fot_mdri2^2)*dm_var_mdri2)   +
      ((fot_frr1^2)*dm_var_frr1)     + ((fot_frr2^2)*dm_var_frr2)
  } else {
    stop("specify BMest")
  }
}

BS_SPREADbyDM <- function (bsdm_spread, bs_count, bsvec, dm_sd)  { # possibly be replace by options from "test_BSbyDM
  spread_seq <- (rnorm(bsdm_spread, mean(bsvec), dm_sd))-mean(bsvec)
  BS_spread <- vector(length=bsdm_spread*bs_count)
  for (i in c(1:bs_count)) {
    BS_spread[(i*bsdm_spread-bsdm_spread+1):(i*bsdm_spread)] <- spread_seq + bsvec[i]
  }
  CIlo <- quantile (BS_spread, 0.025)
  CIup <- quantile (BS_spread, 0.975)
  return(c(CIlo, CIup))
}


#' Incidence and incidence difference from Trinomial Counts / Prevalences of HIV and recency
#'
#' @param BS_Count Add details on each prarameter
#' @return Incidence
#' @return RSE Incidence
#' @return etc.
#' @details
#'
#' general details
#'
#' @examples
#' example
#' @export
recencyI <- function (BS_Count=10000,
                      BSDM_spread=1000,
                      BS_Vars= NULL,
                      BMest="sameTest",
                      PrevH, RSE_PrevH, PrevR, RSE_PrevR,
                      MDRI, RSE_MDRI, FRR, RSE_FRR,
                      BigT=730, Covar_HR=0)

{
  stopifnot (PrevH<=1     & PrevH>=0)
  stopifnot (PrevR<=1     & PrevR>=0)
  stopifnot (RSE_PrevH<=1 & RSE_PrevH>=0)
  stopifnot (RSE_PrevR<=1 & RSE_PrevR>=0)
  stopifnot (MDRI>=0)
  stopifnot (RSE_MDRI<=1  & RSE_MDRI>=0)
  stopifnot (FRR<=1       & FRR>=0)
  stopifnot (RSE_FRR<=1   & RSE_FRR>=0)

  if (length(MDRI)>length(FRR)) {stop("number of inputs for MDRI is larger than number of inputs for MDRI")}
  if (BigT<=182)                {warning ("BigT is smaller than half a year")}

  no_s <- length(PrevH)
  if (length(MDRI)==1)     {MDRI <- rep(MDRI, times=no_s)}         else {MDRI=MDRI}
  if (length(FRR)==1)      {FRR  <- rep(FRR, times=no_s)}          else {FRR=FRR}
  if (length(RSE_MDRI)==1) {RSE_MDRI <- rep(RSE_MDRI, times=no_s)} else {RSE_MDRI=RSE_MDRI}
  if (length(RSE_FRR)==1)  {RSE_FRR  <- rep(RSE_FRR, times=no_s)}  else {RSE_FRR=RSE_FRR}
  if (length(Covar_HR)==1) {Covar_HR <- rep(Covar_HR, times=no_s)} else {Covar_HR=Covar_HR}
  stopifnot(no_s==length(PrevR)   & no_s==length(RSE_PrevH) & no_s==length(PrevR) &
            no_s==length(MDRI)    & no_s==length(RSE_MDRI)  & no_s==length(FRR)   &
            no_s==length(RSE_FRR) & no_s==length(Covar_HR)  & length(BigT)==1)

  MDRI<-MDRI/365.25
  BigT<-BigT/365.25

  I_Est <- I_EST(prevH=PrevH, prevR=PrevR, mdri=MDRI, frr=FRR, bigt=BigT)

  deltaI_Est_Mat<-matrix(ncol=no_s, nrow=no_s)
  for (i in c(1:no_s)) {
    for (j in c(1:no_s)) {
      deltaI_Est_Mat[j,i]<-I_Est[i]-I_Est[j]
    }
  }
  deltaI_Est_Vec <- as.vector(deltaI_Est_Mat)

  if (is.element("PrevH", BS_Vars)) {
    BS_Var_PrevH <- (RSE_PrevH*PrevH)^2
    DM_Var_PrevH <- rep(0, times=no_s)
  } else {
    BS_Var_PrevH <- rep(0, times=no_s)
    DM_Var_PrevH <- (RSE_PrevH*PrevH)^2
  }
  if (is.element("PrevR", BS_Vars)) {
    BS_Var_PrevR <- (RSE_PrevR*PrevR)^2
    DM_Var_PrevR <- rep(0, times=no_s)
  } else {
    BS_Var_PrevR <- rep(0, times=no_s)
    DM_Var_PrevR <- (RSE_PrevR*PrevR)^2
  }
  if (is.element("MDRI", BS_Vars))  {
    BS_Var_MDRI  <- (MDRI*RSE_MDRI)^2
    DM_Var_MDRI  <- rep(0, times=no_s)
  } else {
    BS_Var_MDRI  <- rep(0, times=no_s)
    DM_Var_MDRI  <- (MDRI*RSE_MDRI)^2
  }
  if (is.element("FRR", BS_Vars))   {
    BS_Var_FRR   <- (FRR*RSE_FRR)^2
    DM_Var_FRR   <- rep(0, times=no_s)
  } else {
    BS_Var_FRR   <- rep(0, times=no_s)
    DM_Var_FRR   <- (FRR*RSE_FRR)^2
  }

  if (length(BS_Vars)!=0) {
    I_BSMat <- matrix(nrow=BS_Count, ncol=no_s)
    BS_RootEstMat <- matrix(nrow=BS_Count, ncol=no_s*4)
    BS_Var_I <- vector(length=no_s)
    for (i in c(1:no_s)) {
      BS_RootEstMat [,(i*4-3):(i*4)] <- BS_SURVEY_ESTS  (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i], frr=FRR[i], bs_count=BS_Count,
                                                         bs_var_prevH=BS_Var_PrevH[i], bs_var_prevR=BS_Var_PrevR[i],
                                                         bs_var_mdri=BS_Var_MDRI[i], bs_var_frr=BS_Var_FRR[i],
                                                         covar_HR=Covar_HR[i])
    }
      if ((BMest=="sameTest"| BMest=="FRR.indep")  & (is.element("MDRI", BS_Vars))) {
        BS_RootEstMat[,(i*4-1)] <- BS_RootEstMat[,3]
      }
      if (BMest=="sameTest" & (is.element("FRR", BS_Vars))) {
        BS_RootEstMat[,(i*4)] <- BS_RootEstMat[,4]
      }
    for (i in c(1:no_s)) {
      I_BSVec <- I_EST(prevH=BS_RootEstMat[,(i*4-3)], prevR=BS_RootEstMat[,(i*4-2)],
                       mdri=BS_RootEstMat[,(i*4-1)],  frr=BS_RootEstMat[,(i*4)], bigt=BigT)
      I_BSMat[,i] <-I_BSVec
      BS_Var_I[i] <- var(I_BSMat[,i])
    }

    deltaI_BSMat <- matrix(nrow=BS_Count, ncol=(no_s^2))
    for (i in c(1:no_s)){
      for (j in c(1:no_s)) {
        deltaI_BSMat[,(i*no_s-(no_s-j))] <- I_BSMat[,i]-I_BSMat[,j]
        I_BSMat[,i] <- sort(I_BSMat[,i], decreasing=FALSE)
      }
    }
    BS_Var_deltaI <- vector(length=no_s^2)
    for (i in c(1:(no_s^2))) {
      deltaI_BSMat[,i] <- sort(deltaI_BSMat[,i], decreasing=FALSE)
      BS_Var_deltaI[i] <- var(deltaI_BSMat[,i])
    }
  } else { # if BS_Vars is NULL
    BS_Var_I=rep(0, times=no_s)
    BS_Var_deltaI=rep(0, times=no_s^2)
  }

  if (length(BS_Vars)!=4) {
    fot_Mat  <- matrix (nrow=no_s, ncol=4)
    DM_Var_I <- vector(length=no_s)
    for (i in c(1:no_s)) {
      fot_Mat[i,] <- DM_FirstOrderTerms (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i], frr=FRR[i], bigt=BigT)
      DM_Var_I[i] <- (fot_Mat[i,1]^2)*DM_Var_PrevH[i] + (fot_Mat[i,2]^2)*DM_Var_PrevR[i] +
                     (fot_Mat[i,3]^2)*DM_Var_MDRI[i]  + (fot_Mat[i,4]^2)*DM_Var_FRR[i]
    }

    DM_Var_deltaI <- vector(length=no_s^2)
    for (i in c(1:no_s)) {
      for (j in c(1:no_s)) {
        DM_Var_deltaI[i*no_s-(no_s-j)] <- DM_VAR_deltaI (bmest=BMest, fot_prevH1=fot_Mat[i,1], fot_prevH2=fot_Mat[j,1],
                                                         fot_prevR1=fot_Mat[i,2],   fot_prevR2=fot_Mat[j,2],
                                                         fot_mdri1=fot_Mat[i,3],    fot_mdri2=fot_Mat[j,3],
                                                         fot_frr1=fot_Mat[i,4],     fot_frr2=fot_Mat[j,4],
                                                         dm_var_prevH1=DM_Var_PrevH[i], dm_var_prevH2=DM_Var_PrevH[j],
                                                         dm_var_prevR1=DM_Var_PrevR[i], dm_var_prevR2=DM_Var_PrevR[j],
                                                         dm_var_mdri1=DM_Var_MDRI[i],   dm_var_mdri2=DM_Var_MDRI[j],
                                                         dm_var_frr1=DM_Var_FRR[i],     dm_var_frr2=DM_Var_FRR[j])
      }
    }


  } else { # if BS is performed for all variables
    DM_Var_I <- rep(0, times=no_s)
    DM_Var_deltaI <- rep(0, times=no_s^2)
  }

  DM_SD_I      <- sqrt(DM_Var_I)
  DM_SD_deltaI <- sqrt(DM_Var_deltaI)
  Var_I        <- BS_Var_I + DM_Var_I
  RSE_I        <- sqrt(Var_I)/I_Est
  SD_I         <- sqrt(Var_I)
  Var_deltaI   <- DM_Var_deltaI +  BS_Var_deltaI
  RSE_deltaI   <- sqrt(Var_deltaI)/abs(deltaI_Est_Vec)
  SD_deltaI    <- sqrt(Var_deltaI)

 CI_BSandDM <- function (BSMat, DM_SD, Est) {
    if (sum(DM_Var_I)>0)  {
       if (sum(BS_Var_I)>0)  {
         for (i in c(1:ncol(BSMat))) {
           CI_Mat[i,] <- BS_SPREADbyDM(bsdm_spread=BSDM_spread, bs_count=BS_Count, bsvec=BSMat[,i], dm_sd=DM_SD[i])
         }
       } else {
         for (i in c(1:length(Est))) {
               CI_Mat[i,1] <- qnorm(0.025, mean=Est[i], sd=DM_SD[i])
               CI_Mat[i,2] <- qnorm(0.975, mean=Est[i], sd=DM_SD[i])
         }
       }
    } else {
      for (i in c(1:ncol(BSMat))) {
             CI_Mat[i,1] <- quantile(BSMat[,i], 0.025)
             CI_Mat[i,2] <- quantile(BSMat[,i], 0.975)
      }
    }
    return (CI_Mat)
  }

  CI_Mat <- matrix(nrow=no_s, ncol=2)
  CI_I_Mat <- CI_BSandDM (BSMat=I_BSMat, DM_SD=DM_SD_I, Est=I_Est)
  CI_Mat <- matrix(nrow=no_s^2, ncol=2)
  CI_deltaI_Mat <- CI_BSandDM (BSMat=deltaI_BSMat, DM_SD=DM_SD_deltaI, Est=deltaI_Est_Vec)

  p_value <- pnorm((-abs(deltaI_Est_Vec)/SD_deltaI), mean=0, sd=1)*2

  survey_no <- vector(length=no_s^2)
  out_I_Est <- vector(length=no_s^2)
  out_RSE_I <- vector(length=no_s^2)
  out_CI_I_lo  <- vector(length=no_s^2)
  out_CI_I_up <- vector(length=no_s^2)
  delta_code <- vector(length=no_s^2)
  for (i in c(1:no_s)) {
    survey_no  [(i*no_s-(no_s-1)):(i*no_s)] <- c(i, rep("", times=(no_s-1)))
    out_I_Est  [(i*no_s-(no_s-1)):(i*no_s)] <- c(round(I_Est[i], digit=5), rep("", times=(no_s-1)))
    out_RSE_I  [(i*no_s-(no_s-1)):(i*no_s)] <- c(round(RSE_I[i], digit=5), rep("", times=(no_s-1)))
    out_CI_I_lo[(i*no_s-(no_s-1)):(i*no_s)] <- c(round(CI_I_Mat[i,1], digit=5), rep("", times=(no_s-1)))
    out_CI_I_up[(i*no_s-(no_s-1)):(i*no_s)] <- c(round(CI_I_Mat[i,2], digit=5), rep("", times=(no_s-1)))
    for (j in c(1:no_s)) {
      delta_code [(i*no_s-(no_s-j))] <- paste(i, j, sep=" vs ")
    }
  }

  for (i in c(1:no_s)) {
    deltaI_Est_Vec[(i*no_s-(no_s-i))]<-NA
    RSE_deltaI[(i*no_s-(no_s-i))]<-NA
    p_value[(i*no_s-(no_s-i))]<-NA
    for (j in c(1:2)) {
      CI_deltaI_Mat[(i*no_s-(no_s-i)),j] <- NA
    }
  }

  out_deltaI_Est <- round(deltaI_Est_Vec, digit=5)
  out_RSE_deltaI <- round(RSE_deltaI, digit=5)
  out_p_value <- round(p_value, digit=5)
  out_CI_deltaI_Mat <- round(CI_deltaI_Mat, digit=5)

  if (length(I_Est)==1) {
    output <- data.frame ("Incidence"=out_I_Est,
                          "CI lo"=out_CI_I_lo,
                          "CI up"=out_CI_I_up,
                          "RSE"=out_RSE_I) } else {
    output <- data.frame ("survey"=survey_no,
                          "Incidence"=out_I_Est,
                          "CI lo"=out_CI_I_lo,
                          "CI up"=out_CI_I_up,
                          "RSE"=out_RSE_I,
                          "compare"=delta_code,
                          "Diff"=out_deltaI_Est,
                          "CI Diff lo"=out_CI_deltaI_Mat[,1],
                          "CI Diff up"=out_CI_deltaI_Mat[,2],
                          "RSE Diff"=out_RSE_deltaI,
                          "p-value"=out_p_value) }

  return (output)
}









probs<-prevBYcounts (N=c(5000,5000), N_H=c(1000,1000), N_testR=c(1000,1000), N_R=c(100,70))

################### == Call - DM only ==######################################################
recencyI  (BS_Count=10000,
           BSDM_spread=1000,
           BS_Vars=NULL,
           BMest="sameTest",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)
##############################################################################################

################### == Call - BS only ==######################################################
recencyI  (BS_Count=10000,
           BSDM_spread=1000,
           BS_Vars=c("PrevH", "PrevR", "MDRI","FRR"),
           BMest="sameTest",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)
##############################################################################################

################### == Call - DM & BS combinded ==############################################
recencyI  (BS_Count=10000,
           BSDM_spread=1000,
           BS_Vars=c("MDRI","FRR"),
           BMest="sameTest",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)
##############################################################################################

################### == Only 1 Dataset == Call - DM & BS combinded ==##########################
probs<-prevBYcounts (N=5000, N_H=1000, N_testR=1000, N_R=100)

recencyI  (BS_Count=10000,
           BSDM_spread=1000,
           BS_Vars=c("MDRI","FRR"),
           BMest="sameTest",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)
##############################################################################################

########== Check with spread sheets """"""sameTest""""""==###########
probs<-prevBYcounts (N=c(5000,5000,3000), N_H=c(1000,1000,1000), N_testR=c(1000,1000,900), N_R=c(100,70,120))

recencyI  (BS_Count=10000,
           BSDM_spread=1000,
           BS_Vars=NULL,
           BMest="sameTest",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)

########== Check with spread sheets """"""FRR.indep""""""==###########
probs<-prevBYcounts (N=c(5000,5000,3000), N_H=c(1000,1000,1000), N_testR=c(1000,1000,900), N_R=c(100,70,120))

recencyI  (BS_Count=10000,
           BSDM_spread=1000,
           BS_Vars=NULL,
           BMest="FRR.indep",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)

########== Check with spread sheets """"""MDRI.FRR.indep""""""==###########
probs<-prevBYcounts (N=c(5000,5000,3000), N_H=c(1000,1000,1000), N_testR=c(1000,1000,900), N_R=c(100,70,120))

recencyI  (BS_Count=10000,
           BSDM_spread=1000,
           BS_Vars=NULL,
           BMest="MDRI.FRR.indep",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)

########== bootstrap control by FRR=0 ==###########
probs<-prevBYcounts (N=c(5000,5000,3000), N_H=c(1000,1000,1000), N_testR=c(1000,1000,900), N_R=c(100,70,120))

recencyI  (BS_Count=10000,
           BSDM_spread=1000,
           BS_Vars=c("PrevH","PrevR","MDRI","FRR"),
           BMest="MDRI.FRR.indep",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0, RSE_FRR=0,
           BigT=730)


#   ## control function with I_EST leaving out FRR ####
# BS_Count=10000
# BSDM_spread=1000
# BS_Vars=c("PrevH","PrevR","MDRI","FRR")
# BMest="MDRI.FRR.indep"
# PrevH=probs[,1]
# RSE_PrevH=probs[,3]
# PrevR=probs[,2]
# RSE_PrevR=probs[,4]
# MDRI=c(200/365.25,200/365.25,200/365.25)
# RSE_MDRI=c(0.05,0.05,0.05)
# Covar_HR=c(0,0,0)
# BigT=730/365.25
#
# if 	(is.element("PrevH", BS_Vars)) {BS_Var_PrevH <- (RSE_PrevH*PrevH)^2
# DM_Var_PrevH <- rep(0, times=no_s)}  else
# {BS_Var_PrevH <- rep(0, times=no_s)
# DM_Var_PrevH <- (RSE_PrevH*PrevH)^2}
# if 	(is.element("PrevR", BS_Vars)) {BS_Var_PrevR <- (RSE_PrevR*PrevR)^2
# DM_Var_PrevR <- rep(0, times=no_s)}  else
# {BS_Var_PrevR <- rep(0, times=no_s)
# DM_Var_PrevR <- (RSE_PrevR*PrevR)^2}
# if 	(is.element("MDRI", BS_Vars))  {BS_Var_MDRI  <- (MDRI*RSE_MDRI)^2
# DM_Var_MDRI  <- rep(0, times=no_s)}   else
# {BS_Var_MDRI  <- rep(0, times=no_s)
# DM_Var_MDRI  <- (MDRI*RSE_MDRI)^2}
# if 	(is.element("FRR", BS_Vars))   {BS_Var_FRR   <- (FRR*RSE_FRR)^2
# DM_Var_FRR   <- rep(0, times=no_s)}    else
# {BS_Var_FRR   <- rep(0, times=no_s)
# DM_Var_FRR   <- (FRR*RSE_FRR)^2}
#
# I_EST2 <- function (prevH, prevR, mdri)
# {prevH*prevR/((1 - prevH)*mdri)}
# I_Est2 <- I_EST2(prevH=PrevH, prevR=PrevR, mdri=MDRI)
# I_Est2
#
# BOOTSTRAP <- function(prevH, prevR, mdri, bs_var_prevH, bs_var_prevR, bs_var_mdri, covar_HR)
# {require(MASS)
#   Mu    <- c(prevH, prevR, mdri)
#   sigma <- matrix( c(bs_var_prevH,  covar_HR,     0,
#                      covar_HR,      bs_var_prevR, 0,
#                      0,             0,            bs_var_mdri ),
#                    nrow=3, ncol=3)
#   BS_RootEst <- mvrnorm(n=BS_Count, mu=Mu, Sigma=sigma, empirical=TRUE)
#   return (BS_RootEst)}
#
# no_s <- length(PrevH)
#   I_BSMat <- matrix(nrow=BS_Count, ncol=no_s)
#   BS_RootEstMat <- matrix(nrow=BS_Count, ncol=no_s*3)
#   BS_Var_I <- vector(length=no_s)
#   for (i in c(1:no_s)) {
#     BS_RootEstMat [,(i*3-2):(i*3)] <- BS_SURVEY_ESTS  (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i],
#                                                  bs_var_prevH=BS_Var_PrevH[i], bs_var_prevR=BS_Var_PrevR[i],
#                                                  bs_var_mdri=BS_Var_MDRI[i],
#                                                  covar_HR=Covar_HR[i])
#     if ((BMest=="sameTest"| BMest=="FRR.indep")  & (is.element("MDRI", BS_Vars))) {
#       BS_RootEstMat[,(i*3-1)] <- BS_RootEstMat[,2] }
#     if (BMest=="sameTest" & (is.element("FRR", BS_Vars))) {
#       BS_RootEstMat[,(i*3)] <- BS_RootEstMat[,3] }
#     I_BSVec <- I_EST2(prevH=BS_RootEstMat[,(i*3-2)], prevR=BS_RootEstMat[,(i*3-1)],
#                      mdri=BS_RootEstMat[,(i*3)])
#     I_BSMat[,i] <-I_BSVec
#     BS_Var_I[i] <- var(I_BSMat[,i]) }
#
#   RSE_I2 <- sqrt(BS_Var_I)/I_Est2
#











