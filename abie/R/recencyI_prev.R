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





#' Prevalence and Relative Standard Errors by Counts
#'
#' @param N Counts of total survey sample size(s) (vector/integer).
#' @param N_H Number of HIV positive found in survey(s) (vector/integer).
#' @param N_testR Number tested for recency in survey(s) (vector/integer).
#' @param N_R Number of recent cases in survey(s) (vector/integer).
#' @param DE_H Design effect of HIV prevalence test (vector/numeric), greater than or equal to 1. If multiple surveys are entered but only one design effect is specified, function assumes entered design effect is identical for both surveys.
#' @param DE_R Design effect of recency test (vector/numeric), greater than or equal to 1. If multiple surveys are entered but only one design effect is specified, function assumes entered design effect is identical for both surveys.
#' @return Prevalences and relative standard errors. Design effects are assumed negligible unless user specifies otherwise.
#'
#' @examples
#' prevBYcounts(N = 5000, N_H = 1000, N_testR = 1000, N_R = 70, DE_R = 1.1)
#'
#' prevBYcounts (N = c(5000,5000), N_H = c(1000,1000), N_testR = c(1000,1000),
#' N_R = c(100,70), DE_H = 1, DE_R = c(1,1.1))
#' @export
#'
prevBYcounts <- function(N, N_H, N_testR, N_R, DE_H=1, DE_R=1) {
  if (sum(N_H>N)>0 | sum(N_R>N)>0) {
    stop("sample subset larger than total sample")
  }
  if (sum(N_testR<N_R)>0) {
    stop("counts of recency tested less than counts of those found to be recently infected")
  }
  if (sum(DE_H<1) >0 | sum(DE_R<1) >0) {
    stop("Design effects must be >=1")
  }
  if(sum(N_testR/N_H<.95)>0){
    warning("Less than 95% of HIV-positive subjects have recency test results available")
  }

  no_s <- length(N)
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


#SIMPLE FORMULA FOR INCIDENCE
I_EST <- function (prevH, prevR, mdri, frr, bigt){
                   prevH*(prevR - frr)/((1 - prevH)*(mdri-frr*bigt))
}


#BOOTSTRAPPING EDF OF INPUT PARAMETERS TO RECENCYI() FUNCTION
BS_SURVEY_ESTS <- function (prevH, prevR, mdri, frr, bs_count,
                            bs_var_prevH, bs_var_prevR, bs_var_mdri, bs_var_frr, covar_HR) {
  Mu    <- c(prevH, prevR, mdri, frr)
  sigma <- matrix( c(bs_var_prevH,  covar_HR,     0,           0,
                     covar_HR,      bs_var_prevR, 0,           0,
                     0,             0,            bs_var_mdri, 0,
                     0,             0,            0,           bs_var_frr),
                 nrow=4, ncol=4)
#bs_var_prevH,   bs_var_prevR, and so on are...?empirical, observed variance of variable?
#returns bootstraps of prevH, prevR, mdri, frr
  BS_RootEst <- MASS::mvrnorm(n=bs_count, mu=Mu, Sigma=sigma, empirical=TRUE)
  return (BS_RootEst)
}


#INPUT IS PROPORTIONS, TEST CHARACTERISTICS, OUTPUT IS PARTIAL DERIVATIVES WITH RESPECT TO (WRT) EACH INPUT VARIABLE
DM_FirstOrderTerms <- function (prevH, prevR, mdri, frr, bigt)   {
  fot_prevH <- (prevR-frr)/(((1-prevH)^2)*(mdri-frr*bigt)) #E.G. d(I)/d(P_H)
  fot_prevR <- prevH/((1-prevH)*(mdri-frr*bigt))
  fot_mdri  <- (frr*prevH-prevR*prevH)/((1-prevH)*((mdri-frr*bigt)^2))
  fot_frr   <- (prevH*(bigt*prevR-mdri))/((1-prevH)*((mdri-frr*bigt)^2))
  return (c(fot_prevH, fot_prevR, fot_mdri, fot_frr))
}


#TAKES TESTING SCHEME (SAME ASSAY, DIFFERENT, ETC.) AND RETURNS VAR[diff in I].
DM_VAR_deltaI <- function (BMest, fot_prevH1, fot_prevR1, fot_mdri1, fot_frr1,
                           fot_prevH2, fot_prevR2, fot_mdri2, fot_frr2,
                           dm_var_prevH1, dm_var_prevR1, dm_var_mdri1, dm_var_frr1,
                           dm_var_prevH2, dm_var_prevR2, dm_var_mdri2, dm_var_frr2) {
  if(BMest=="same.test"){
    DM_Var_deltaI <- ((fot_prevH1^2)*dm_var_prevH1) + ((fot_prevH2^2)*dm_var_prevH2) +
      ((fot_prevR1^2)*dm_var_prevR1) + ((fot_prevR2^2)*dm_var_prevR2) +
      ((fot_mdri1-fot_mdri2)^2*dm_var_mdri1) +
      ((fot_frr1-fot_frr2)^2*dm_var_frr1)
  }
  else if(BMest=="FRR.indep") {
    DM_Var_deltaI <- ((fot_prevH1^2)*dm_var_prevH1) + ((fot_prevH2^2)*dm_var_prevH2) +
      ((fot_prevR1^2)*dm_var_prevR1) + ((fot_prevR2^2)*dm_var_prevR2) +
      ((fot_mdri1-fot_mdri2)^2*dm_var_mdri1) +
      ((fot_frr1^2)*dm_var_frr1)     + ((fot_frr2^2)*dm_var_frr2)
  }
  else if(BMest=="MDRI.FRR.indep") {
    DM_Var_deltaI <- ((fot_prevH1^2)*dm_var_prevH1) + ((fot_prevH2^2)*dm_var_prevH2) +
      ((fot_prevR1^2)*dm_var_prevR1) + ((fot_prevR2^2)*dm_var_prevR2) +
      ((fot_mdri1^2)*dm_var_mdri1)   + ((fot_mdri2^2)*dm_var_mdri2)   +
      ((fot_frr1^2)*dm_var_frr1)     + ((fot_frr2^2)*dm_var_frr2)
  } else {
    stop("specify BMest")
  }
}



#TAKES TESTING SCHEME (SAME ASSAY, DIFFERENT, ETC.) AND RETURNS VAR[diff in I] at infinite sample size.
DM_VAR_deltaI.infSS<-function(BMest,
                              fot_mdri1, fot_frr1,
                              fot_mdri2, fot_frr2,
                              dm_var_mdri1, dm_var_frr1,
                              dm_var_mdri2, dm_var_frr2){
  if(BMest=="same.test"){
    DM_Var_deltaI <-
      ((fot_mdri1-fot_mdri2)^2*dm_var_mdri1) +
      ((fot_frr1-fot_frr2)^2*dm_var_frr1)
  }
  else if(BMest=="FRR.indep") {
    DM_Var_deltaI <-
      ((fot_mdri1-fot_mdri2)^2*dm_var_mdri1) +
      ((fot_frr1^2)*dm_var_frr1)     + ((fot_frr2^2)*dm_var_frr2)
  }
  else if(BMest=="MDRI.FRR.indep") {
    DM_Var_deltaI <-
      ((fot_mdri1^2)*dm_var_mdri1)   + ((fot_mdri2^2)*dm_var_mdri2)   +
      ((fot_frr1^2)*dm_var_frr1)     + ((fot_frr2^2)*dm_var_frr2)
  } else {
    stop("specify BMest")
  }
}




#example to get function working:

# probs<-prevBYcounts (N=c(5000,5000), N_H=c(1000,1000), N_testR=c(1000,1000), N_R=c(100,70))
# probs
#
# PrevH = probs[,1]
# RSE_PrevH = probs[,3]
# PrevR = probs[,2]
# RSE_PrevR = probs[,4]
#
#
#
# BS_Count = 10000
# Boot = FALSE
# BMest = "MDRI.FRR.indep"
# MDRI = 200
# RSE_MDRI = 0.05
# FRR = c(0.01,0.009)
# RSE_FRR = 0.2
# BigT = 730
# Covar_HR=0
# alpha=0.05
#
#
# recencyI(PrevH=probs[,1], RSE_PrevH=probs[,3], PrevR=probs[,2], RSE_PrevR=probs[,4], Boot = FALSE, BS_Count = 10000, alpha = 0.05,
# BMest = "MDRI.FRR.indep", MDRI=MDRI, RSE_MDRI=RSE_MDRI, FRR=FRR, RSE_FRR=RSE_FRR, BigT = c(730,720), Covar_HR = 0)
#


#' Incidence and incidence difference statistics from trinomial prevalences of HIV and recency
#'
#' @param PrevH Prevalence of HIV (vector/integer).
#' @param RSE_PrevH Relative Standard Error (RSE) of estimate for population prevalence of HIV (vector/integer).
#' @param PrevR Proportion of persons found to be 'recent' by biomarker assay among total persons found positive for HIV (vector/integer).
#' @param RSE_PrevR Relative Standard Error (RSE) of estimate for population proportion of those testing positive for HIV who have been infected recently (vector/integer).
#' @param BS_Count Specifies number of bootstrap samples for bootstrapped confidence intervals of incidence.
#' @param Boot True/False variable indicating whether variance of point estimates is to be calculated by Empirical Bootstrapping (TRUE) or Delta Method (FALSE), the default setting.
#' @param BMest Biomarker estimation by one the 3 options "same.test"(=default), "FRR.indep", "MDRI.FRR.indep" (string).
#' @param MDRI mean duration of recent infection [days] (vector/integer).
#' @param RSE_MDRI Relative standard error of MDRI [days] (vector/integer).
#' @param FRR False recent rate (vector/integer).
#' @param RSE_FRR Relative standard error of FRR (vector/integer).
#' @param BigT post-infection time cut-off true vs false recent [days] default 730 days (integer).
#' @param Covar_HR Covariance of probability of being positive and being categorized recent from survey (vector/integer). Note that as the variances of PrevH and PrevR are often quite small, only a suitably commensurate covariance will enable the inversion of the bootstrap covariance matrix for random number generation to proceed without error.
#' @return Incidence estimate, confidence interval, relative standard error of estimate and of assay characteristics MDRI and FRR. If multiple surveys are entered, function returns said results, as well as estimates of incidence  differences, confidence intervals of differences, difference relative standard errors, and p-values testing the hypothesis that the difference in incidence measures are zero. Theoretical relative standard error of incidence and incidence difference at infinite sample size is returned only if Boot=FALSE, as that calculation relies on the asymptotic behavior of components of the Delta method approximation, and is not calculable from bootstrapped values.
#' @details Implements assay-based incidence estimation through cross-sectional prevalence and recency of infection tests as described by Kassanjee, et al. "A new general biomarker-based incidence estimator," \emph{Epidemiology} (2012). Function parameters must be specified to include assay test characteristics and survey results as proportions. Confidence intervals are computed via Delta method approximation, except when Boot=TRUE is specified, in which case confidence intervals are generated by empirical bootstrap resampling. Inputs must be in appropriate ranges for appropriate units. Extreme input values may make calculation impossible, and if entered will elicit error notices.
#'
#' @examples
#' recencyI(PrevH = 0.20, RSE_PrevH = 0.028, PrevR = 0.10, RSE_PrevR = 0.094,
#' BS_Count = 10000, Boot = TRUE, BMest = "same.test", MDRI = 200, RSE_MDRI = 0.05,
#' FRR = 0.01, RSE_FRR = 0.2, BigT = 730)
#'
#'
#'recencyI(PrevH = c(0.20,0.21,0.18), RSE_PrevH = c(0.028,0.03,0.022),
#'PrevR = c(0.10,0.13,0.12), RSE_PrevR = c(0.094,0.095,0.05),
#'BS_Count = 10000, Boot = FALSE, BMest = "MDRI.FRR.indep", MDRI = 200,
#'RSE_MDRI = 0.05, FRR = c(0.01,0.009,0.02), RSE_FRR = 0.2, BigT = 730)
#' @export

recencyI <- function (PrevH, RSE_PrevH, PrevR, RSE_PrevR,
                      Boot = FALSE, BS_Count = 10000, alpha = 0.05,
                      BMest = "same.test", MDRI, RSE_MDRI, FRR, RSE_FRR,
                      BigT = 730, Covar_HR = 0){

  stopifnot (PrevH<=1     & PrevH>=0)
  stopifnot (PrevR<=1     & PrevR>=0)
  stopifnot (RSE_PrevH<=1 & RSE_PrevH>=0)
  stopifnot (RSE_PrevR<=1 & RSE_PrevR>=0)
  stopifnot (MDRI>=0)
  stopifnot (RSE_MDRI<=1  & RSE_MDRI>=0)
  stopifnot (FRR<=1       & FRR>=0)
  stopifnot (RSE_FRR<=1   & RSE_FRR>=0)

  if(sum(BMest==c("same.test", "FRR.indep", "MDRI.FRR.indep"))==0){
    stop("BMest option must be same.test, FRR.indep, or MDRI.FRR.indep")
  }

  if(BS_Count<=0 & Boot==TRUE){
    stop("Bootstrap samples count must be positive integer")
  }

  if (BMest!="MDRI.FRR.indep" & length(MDRI)>length(FRR)) {stop("number of inputs for MDRI is larger than number of inputs for FRR")}

  if(sum(MDRI<90)>0){
    warning("Estimated MDRI less than 90 days")
  }

  if(sum(FRR>0.10)>0){
    warning("Estimated FRR is greater than 10%")
  }
  if(sum(FRR==0)>0){
    warning("Zero estimated FRR")
  }
  if(sum(RSE_FRR>0.30)>0){
    warning("RSE of estimated FRR is greater than 30%")
  }
  if(sum(RSE_FRR<0.05)>0){
    warning("RSE of estimated FRR is less than 5%")
  }
  if(sum(RSE_MDRI<0.01)>0){
    warning("RSE of estimated MDRI is less than 1%")
  }
  if(sum(RSE_MDRI>0.30)>0){
    warning("RSE of estimated MDRI is greater than 30%")
  }

  if (sum(MDRI>BigT)>0) {stop("MDRI cannot be greater than BigT")}

  if (sum(BigT<=182)>0) {warning ("BigT is smaller than half a year")}

  no_s <- length(PrevH) #dimension of inputs (number of surveys)

  if (length(MDRI)==1)     {MDRI <- rep(MDRI, times=no_s)}         else {MDRI=MDRI}
  if (length(FRR)==1)      {FRR  <- rep(FRR, times=no_s)}          else {FRR=FRR}
  if (length(RSE_MDRI)==1) {RSE_MDRI <- rep(RSE_MDRI, times=no_s)} else {RSE_MDRI=RSE_MDRI}
  if (length(RSE_FRR)==1)  {RSE_FRR  <- rep(RSE_FRR, times=no_s)}  else {RSE_FRR=RSE_FRR}
  if (length(Covar_HR)==1) {Covar_HR <- rep(Covar_HR, times=no_s)} else {Covar_HR=Covar_HR}
  if (length(BigT)>1 & BMest!="MDRI.FRR.indep") stop("More than one BigT specified when only one MDRI")
  if (length(BigT)==1) {BigT <- rep(BigT, times=no_s)} else {BigT=BigT}
  stopifnot(no_s==length(PrevR)   & no_s==length(RSE_PrevH) & no_s==length(PrevR) &
            no_s==length(MDRI)    & no_s==length(RSE_MDRI)  & no_s==length(FRR)   &
            no_s==length(RSE_FRR) & no_s==length(Covar_HR))

  MDRI<-MDRI/365.25
  BigT<-BigT/365.25





  I_Est <- I_EST(prevH=PrevH, prevR=PrevR, mdri=MDRI, frr=FRR, bigt=BigT)

  deltaI_Est_Mat<-matrix(ncol=no_s, nrow=no_s) #empty matrix of dim (# surveys*# surveys)
  for (i in c(1:no_s)) {
    for (j in c(1:no_s)) {
      deltaI_Est_Mat[j,i]<-I_Est[i]-I_Est[j]
    }
  }
#now each element in deltaI_Est_Mat is the difference in I for each survey
#the ij element is I[i]-I[j].


  deltaI_Est_Vec <- as.vector(deltaI_Est_Mat)
#makes above matrix into a vector

#this section defines the empirical variance of each object to be bootsrapped for the mvtrnorm function
  if (Boot==TRUE) {
    BS_Var_PrevH <- (RSE_PrevH*PrevH)^2
    BS_Var_PrevR <- (RSE_PrevR*PrevR)^2
    BS_Var_MDRI  <- (MDRI*RSE_MDRI)^2
    BS_Var_FRR   <- (FRR*RSE_FRR)^2

    DM_Var_PrevH <- rep(0, times=no_s)
    DM_Var_PrevR <- rep(0, times=no_s)
    DM_Var_MDRI  <- rep(0, times=no_s)
    DM_Var_FRR   <- rep(0, times=no_s)
  } else {
    BS_Var_PrevH <- rep(0, times=no_s)
    BS_Var_PrevR <- rep(0, times=no_s)
    BS_Var_MDRI  <- rep(0, times=no_s)
    BS_Var_FRR   <- rep(0, times=no_s)

    DM_Var_PrevH <- (RSE_PrevH*PrevH)^2
    DM_Var_PrevR <- (RSE_PrevR*PrevR)^2
    DM_Var_MDRI  <- (MDRI*RSE_MDRI)^2
    DM_Var_FRR   <- (FRR*RSE_FRR)^2
  }


  if (Boot==TRUE) {
    DM_Var_I <- rep(0, times=no_s)
    DM_Var_deltaI <- rep(0, times=no_s^2)

    I_BSMat <- matrix(nrow=BS_Count, ncol=no_s) #creates empty matrix of dim (#BS samples*#surveys)
    BS_RootEstMat <- matrix(nrow=BS_Count, ncol=no_s*4) #creates empty matrix of dim (#BS samples-by-#surveys*4)
    BS_Var_I <- vector(length=no_s) #vector of length # of surveys

#loop that, for each survey
      for (i in 1:no_s) {
      BS_RootEstMat [,(i*4-3):(i*4)] <- BS_SURVEY_ESTS  (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i], frr=FRR[i], bs_count=BS_Count,
                                                         bs_var_prevH=BS_Var_PrevH[i], bs_var_prevR=BS_Var_PrevR[i],
                                                         bs_var_mdri=BS_Var_MDRI[i], bs_var_frr=BS_Var_FRR[i],
                                                         covar_HR=Covar_HR[i])
    }
    #returns bootstraps of prevH, prevR, mdri, frr as columns, one for each survey, so if survey#=3, then
    #first four columns are for prevH, prevR, mdri, frr as columns, then next four are for those varaibles from second survey, and so on...
    # dim(BS_RootEstMat)

#I THINK THIS BELOW SECTION SHOULD HAVE AN INDEX LOOP RIGHT?? AS IT STANDS THERE'S NOTHING...
    for(i in 1:no_s){
      if ((BMest=="same.test"| BMest=="FRR.indep")  ) {
        BS_RootEstMat[,(i*4-1)] <- BS_RootEstMat[,3]
      }#above is saying, if the test is the same, or FRR.indep (meaning mdri still same) then each column in BS matrix
       #corresponding to mdri will equal the first BS sample of that variable
      if (BMest=="same.test" ) {
        BS_RootEstMat[,(i*4)] <- BS_RootEstMat[,4]
      } #similarly if it's same test then FRR will be recorded as same in there.
}


    for (i in 1:no_s) {
      I_BSVec <- I_EST(prevH=BS_RootEstMat[,(i*4-3)], prevR=BS_RootEstMat[,(i*4-2)],
                       mdri=BS_RootEstMat[,(i*4-1)],  frr=BS_RootEstMat[,(i*4)], bigt=BigT[i])
      I_BSMat[,i] <-I_BSVec
      BS_Var_I[i] <- var(I_BSMat[,i])
    }
#now compute BS I and variance of I from BS samples. dim(I_BSMat) (same rows as number of surveys)


#matrix of bootstrapped differences between surveys
  deltaI_BSMat <- matrix(nrow=BS_Count, ncol=(no_s^2))#empty matrix of dim (BS samples by number surveys squared)
    for (i in c(1:no_s)){
      for (j in c(1:no_s)) {
        deltaI_BSMat[,(i*no_s-(no_s-j))] <- I_BSMat[,i]-I_BSMat[,j]
        I_BSMat[,i] <- sort(I_BSMat[,i], decreasing=FALSE)
      }
    }
#above makes empty matrix of nrow=#BS samples, ncol=#surveys squared, and inputs to each the difference for each BS estimate of incidence

#makes vector of variances of difference in BS estimates of I from different surveys
    BS_Var_deltaI <- vector(length=no_s^2)
    for (i in c(1:(no_s^2))) {
      deltaI_BSMat[,i] <- sort(deltaI_BSMat[,i], decreasing=FALSE)
      BS_Var_deltaI[i] <- var(deltaI_BSMat[,i])
    } #above sorts bootstrap differences and gets an estimate of the variance of each difference (some redundancy here)

#tracks back to line: 'if (Boot==TRUE) {'

  }
  else { # if Boot is FALSE
    BS_Var_I=rep(0, times=no_s) #so if no BS-ing, make the BS variables null
    BS_Var_deltaI=rep(0, times=no_s^2)

    #next few lines make delta method matrix
    fot_Mat  <- matrix (nrow=no_s, ncol=4)
    DM_Var_I <- vector(length=no_s)
    DM_Var_I.infSS <- vector(length=no_s) #ADDED 4.19
    for (i in c(1:no_s)) {
      fot_Mat[i,] <- DM_FirstOrderTerms (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i], frr=FRR[i], bigt=BigT[i])
      #DM_FirstOrderTerms gives fot_prevH, fot_prevR, fot_mdri, fot_frr, so FOT for each prev. survey input
      #and each row of fot_Mat has the FOT for prev. survey i.
      DM_Var_I[i] <- (fot_Mat[i,1]^2)*DM_Var_PrevH[i] + (fot_Mat[i,2]^2)*DM_Var_PrevR[i] +
                     (fot_Mat[i,3]^2)*DM_Var_MDRI[i]  + (fot_Mat[i,4]^2)*DM_Var_FRR[i]
      #DM_Var_I uses the fot's and variances to get the variance for each survey
      DM_Var_I.infSS[i]<-(fot_Mat[i,3]^2)*DM_Var_MDRI[i]  + (fot_Mat[i,4]^2)*DM_Var_FRR[i]
    }

    DM_Var_deltaI <- vector(length=no_s^2) #creates vector of length number of surveys squared.
    DM_Var_deltaI.infSS <- vector(length=no_s^2) #creates vector of length number of surveys squared.
    for (i in c(1:no_s)) {
      for (j in c(1:no_s)) { #starts at 2-1 = 1, goes to 4-0 = 4.
        DM_Var_deltaI[i*no_s-(no_s-j)] <- DM_VAR_deltaI(BMest=BMest, fot_prevH1=fot_Mat[i,1], fot_prevH2=fot_Mat[j,1],
                                                         fot_prevR1=fot_Mat[i,2],   fot_prevR2=fot_Mat[j,2],
                                                         fot_mdri1=fot_Mat[i,3],    fot_mdri2=fot_Mat[j,3],
                                                         fot_frr1=fot_Mat[i,4],     fot_frr2=fot_Mat[j,4],
                                                         dm_var_prevH1=DM_Var_PrevH[i], dm_var_prevH2=DM_Var_PrevH[j],
                                                         dm_var_prevR1=DM_Var_PrevR[i], dm_var_prevR2=DM_Var_PrevR[j],
                                                         dm_var_mdri1=DM_Var_MDRI[i],   dm_var_mdri2=DM_Var_MDRI[j],
                                                         dm_var_frr1=DM_Var_FRR[i],     dm_var_frr2=DM_Var_FRR[j])

        DM_Var_deltaI.infSS[i*no_s-(no_s-j)] <-DM_VAR_deltaI.infSS(BMest=BMest,fot_mdri1=fot_Mat[i,3], fot_frr1=fot_Mat[i,4],
                                                                   fot_mdri2=fot_Mat[j,3], fot_frr2=fot_Mat[j,4],
                                                                   dm_var_mdri1=DM_Var_MDRI[i], dm_var_frr1=DM_Var_FRR[i],
                                                                   dm_var_mdri2=DM_Var_MDRI[j], dm_var_frr2=DM_Var_FRR[j])
        }
    } #creates a vector of length number of surveys squared, each term is delta method variance applied to
    #survey 1 & 1, 1&2, 2&1, and 2&2.
    RSE_I.infSS  <- sqrt(DM_Var_I.infSS)/I_Est #only given if boot=F
    RSE.deltaI.infSS<- sqrt(DM_Var_deltaI.infSS)/abs(deltaI_Est_Vec)
  }


  DM_SD_I      <- sqrt(DM_Var_I)
  DM_SD_deltaI <- sqrt(DM_Var_deltaI)
  Var_I        <- BS_Var_I + DM_Var_I
  RSE_I        <- sqrt(Var_I)/I_Est
  Var_MDRI <- BS_Var_MDRI + DM_Var_MDRI
  Var_FRR  <- BS_Var_FRR + DM_Var_FRR

  # RSE_I.infSS  <- sqrt(DM_Var_I.infSS)/I_Est #only given if boot=F
  # RSE.deltaI.infSS<- sqrt(DM_Var_deltaI.infSS)/abs(deltaI_Est_Vec)
  #put these in the boot=F section
  SD_I         <- sqrt(Var_I)
  Var_deltaI   <- DM_Var_deltaI +  BS_Var_deltaI
  RSE_deltaI   <- sqrt(Var_deltaI)/abs(deltaI_Est_Vec)
  SD_deltaI    <- sqrt(Var_deltaI)

#this function takes bootstrap matrix, SD via delta method, and I estimates, and returns the spread version of those
#IM REWRITING THIS SO THAT IT ONLY TAKES FULL BS OR DELTA-METHOD. Now the method takes estimates and SE and gives CIs
 CI_BSandDM <- function (BSMat, DM_SD, Est) {
    if (sum(DM_Var_I)>0)  {
         for (i in c(1:length(Est))) {
               CI_Mat[i,1] <- qnorm(alpha/2, mean=Est[i], sd=DM_SD[i])
               CI_Mat[i,2] <- qnorm(1-alpha/2, mean=Est[i], sd=DM_SD[i])
         }
       }
    else {
      for (i in c(1:ncol(BSMat))) {
             CI_Mat[i,1] <- quantile(BSMat[,i], alpha/2)
             CI_Mat[i,2] <- quantile(BSMat[,i], 1-alpha/2)
      }
    }
    return (CI_Mat)
  }

  CI_Mat <- matrix(nrow=no_s, ncol=2)
  CI_I_Mat <- CI_BSandDM (BSMat=I_BSMat, DM_SD=DM_SD_I, Est=I_Est) #CIs for incidence for each survey
  CI_Mat <- matrix(nrow=no_s^2, ncol=2)
  CI_deltaI_Mat <- CI_BSandDM (BSMat=deltaI_BSMat, DM_SD=DM_SD_deltaI, Est=deltaI_Est_Vec)
  #above gives CIs for each element in deltaI_Est_Vec, which is a four-tuple vector for 2 surveys

  p_value <- pnorm((-abs(deltaI_Est_Vec)/SD_deltaI), mean=0, sd=1)*2
  if(Boot==F){p_value.infSS <- pnorm((-abs(deltaI_Est_Vec)/sqrt(DM_Var_deltaI.infSS)), mean=0, sd=1)*2}

  #gives 4 p-values for 2 surveys

  survey_no <- vector(length=no_s^2) #gives 1 for 1 survey, 4 for two surveys, etc.
  out_I_Est <- vector(length=no_s^2)
  out_RSE_I <- vector(length=no_s^2)
  out_CI_I_lo  <- vector(length=no_s^2)
  out_CI_I_up <- vector(length=no_s^2)
  delta_code <- vector(length=no_s^2)
  #need to add vector out_RSE_I_infSS here
  if(Boot==F){
  out_RSE_I.infSS <- vector(length=no_s^2) #ADDED 4.19
  out_RSE.deltaI.infSS<- RSE.deltaI.infSS
  }

  for (i in c(1:no_s)) {
    survey_no  [(i*no_s-(no_s-1)):(i*no_s)] <- c(i, rep("", times=(no_s-1)))
    out_I_Est  [(i*no_s-(no_s-1)):(i*no_s)] <- c(round(I_Est[i], digit=5), rep("", times=(no_s-1)))
    out_RSE_I  [(i*no_s-(no_s-1)):(i*no_s)] <- c(round(RSE_I[i], digit=5), rep("", times=(no_s-1)))
    #need to add vector out_RSE_I_infSS here
    if(Boot==F){
    out_RSE_I.infSS  [(i*no_s-(no_s-1)):(i*no_s)] <- c(round(RSE_I.infSS[i], digit=5), rep("", times=(no_s-1)))
    }

    out_CI_I_lo[(i*no_s-(no_s-1)):(i*no_s)] <- c(round(CI_I_Mat[i,1], digit=5), rep("", times=(no_s-1)))
    out_CI_I_up[(i*no_s-(no_s-1)):(i*no_s)] <- c(round(CI_I_Mat[i,2], digit=5), rep("", times=(no_s-1)))
    for (j in c(1:no_s)) {
      delta_code [(i*no_s-(no_s-j))] <- paste(i, j, sep=" vs ")
    }
    }
  #above makes components of output matrix, CIs, empty spaces, etc.

  if(sum(out_RSE_I>0.25)){warning("RSE of incidence estimator greater than 25%")}

  #puts empty spaces in each place of vectors needed (the ends)
  for (i in c(1:no_s)) {
    deltaI_Est_Vec[(i*no_s-(no_s-i))]<-NA

    if(Boot==F){RSE_deltaI[(i*no_s-(no_s-i))]<-NA
    out_RSE.deltaI.infSS[(i*no_s-(no_s-i))]<-NA}
    p_value[(i*no_s-(no_s-i))]<-NA
    if(Boot==F){p_value.infSS[(i*no_s-(no_s-i))]<-NA}
    for (j in c(1:2)) {
      CI_deltaI_Mat[(i*no_s-(no_s-i)),j] <- NA
    }
  }



    for(i in 1:length(out_CI_I_lo)){
      if((out_CI_I_lo[i]!="" & out_CI_I_lo[i]<0) | ( (out_CI_I_lo[i]!="" & out_CI_I_lo[i]>1)))
      {warning("CI out of [0,1] bounds"); break}
    }
  for(i in 1:length(out_CI_I_up)){
    if((out_CI_I_up[i]!="" & out_CI_I_up[i]<0) | ( (out_CI_I_up[i]!="" & out_CI_I_up[i]>1)))
    {warning("CI out of [0,1] bounds"); break}
  }


  out_deltaI_Est <- round(deltaI_Est_Vec, digit=5)

  out_RSE_deltaI <- round(RSE_deltaI, digit=5)

  if(Boot==F){
    out_p_value.infSS <- round(p_value.infSS,5)
    out_RSE.deltaI.infSS <- round(out_RSE.deltaI.infSS, digit=5)
    out_p_value.infSS<-ifelse(out_p_value.infSS<0.001,"<0.0001",out_p_value.infSS)
  }
  out_p_value <- round(p_value,5)
  out_CI_deltaI_Mat <- round(CI_deltaI_Mat, digit=5)

  out_p_value<-ifelse(out_p_value<0.001,"<0.0001",out_p_value)

  MDRI.CI <- round(365.25*data.frame(CI.low=qnorm(alpha/2, mean=MDRI, sd=sqrt(Var_MDRI)), CI.up=qnorm(1-alpha/2, mean=MDRI, sd=sqrt(Var_MDRI))),3)
  FRR.CI <- round(data.frame(CI.low=qnorm(alpha/2, mean=FRR, sd=sqrt(Var_FRR)),CI.up=qnorm(1-alpha/2, mean=FRR, sd=sqrt(Var_FRR))),4)






  if (length(I_Est)==1) {
    output <- list(Incidence.Statistics=data.frame("Incidence"=out_I_Est,
                          "CI low"=out_CI_I_lo,
                          "CI up"=out_CI_I_up,
                          "RSE"=out_RSE_I,
                          "RSE.Inf.SS"=out_RSE_I.infSS), MDRI.CI=MDRI.CI, FRR.CI=FRR.CI)
    } else if(Boot==F){
      Incidence.Statistics=data.frame ("survey"=survey_no,
                                       "Incidence"=out_I_Est,
                                       "CI low"=out_CI_I_lo,
                                       "CI up"=out_CI_I_up,
                                       "RSE"=out_RSE_I,
                                       "RSE.Inf.SS"=out_RSE_I.infSS)
      Incidence.Difference.Statistics=data.frame(
        "compare"=delta_code,
        "Diff"=out_deltaI_Est,
        "CI Diff low"=out_CI_deltaI_Mat[,1],
        "CI Diff up"=out_CI_deltaI_Mat[,2],
        "RSE Diff"=out_RSE_deltaI,
        "RSE Diff Inf.SS"=out_RSE.deltaI.infSS,
        "p-value"=out_p_value,
        "p-value.Inf.SS"=out_p_value.infSS)
      Incidence.Statistics=Incidence.Statistics[which(Incidence.Statistics[,1]!=""),]
      Incidence.Difference.Statistics=Incidence.Difference.Statistics[which(!is.na(Incidence.Difference.Statistics[,2])),]
      row.names(Incidence.Statistics)<-1:nrow(Incidence.Statistics)
      row.names(Incidence.Difference.Statistics)<-1:nrow(Incidence.Difference.Statistics)

    output <- list(Incidence.Statistics=Incidence.Statistics,
                   Incidence.Difference.Statistics=Incidence.Difference.Statistics, MDRI.CI=MDRI.CI, FRR.CI=FRR.CI)
    }else { Incidence.Difference.Statistics = data.frame("compare"=delta_code,
                                                   "Diff"=out_deltaI_Est,
                                                   "CI Diff low"=out_CI_deltaI_Mat[,1],
                                                   "CI Diff up"=out_CI_deltaI_Mat[,2],
                                                   "RSE Diff"=out_RSE_deltaI,
                                                   "p-value"=out_p_value)
      Incidence.Statistics = data.frame("survey"=survey_no,
                                        "Incidence"=out_I_Est,
                                        "CI low"=out_CI_I_lo,
                                        "CI up"=out_CI_I_up,
                                        "RSE"=out_RSE_I)
      Incidence.Statistics = Incidence.Statistics[which(Incidence.Statistics[,1]!=""),]
      Incidence.Difference.Statistics =Incidence.Difference.Statistics[which(!is.na(Incidence.Difference.Statistics[,2])),]
      row.names(Incidence.Statistics)<-1:nrow(Incidence.Statistics)
      row.names(Incidence.Difference.Statistics)<-1:nrow(Incidence.Difference.Statistics)
      output <- list(Incidence.Statistics = Incidence.Statistics,
                     Incidence.Difference.Statistics =Incidence.Difference.Statistics, MDRI.CI=MDRI.CI, FRR.CI=FRR.CI)
    }


  return (output)
}











#' Incidence and incidence difference statistics from trinomial counts of HIV and recency
#'
#' @param N Counts of total survey sample size(s) (vector/integer).
#' @param N_H Number of HIV positive found in survey(s) (vector/integer).
#' @param N_testR Number tested for recency in survey(s) (vector/integer).
#' @param N_R Number of recent cases in survey(s) (vector/integer).
#' @param DE_H Design effect of HIV prevalence test (vector/numeric), greater than or equal to 1. If multiple surveys are entered but only one design effect is specified, function assumes entered design effect is identical for both surveys.
#' @param DE_R Design effect of recency test (vector/numeric), greater than or equal to 1. If multiple surveys are entered but only one design effect is specified, function assumes entered design effect is identical for both surveys.
#' @param BS_Count Specifies number of bootstrap samples for bootstrapped confidence intervals of incidence.
#' @param Boot True/False variable indicating whether variance of point estimates is to be calculated by Empirical Bootstrapping (TRUE) or Delta Method (FALSE), the default setting.
#' @param BMest Biomarker estimation by one the 3 options "same.test"(=default), "FRR.indep", "MDRI.FRR.indep" (string).
#' @param MDRI mean duration of recent infection [days] (vector/integer).
#' @param RSE_MDRI Relative standard error of MDRI [days] (vector/integer).
#' @param FRR False recent rate (vector/integer).
#' @param RSE_FRR Relative standard error of FRR (vector/integer).
#' @param BigT Cut point in days of recency used in biomarker assay to test recency in a given prevalence survey.
#' @param Covar_HR Covariance of probability of being postive and being categorized recent from survey (or as a vector for multiple surveys).
#' @details Implements assay-based incidence estimation through cross-sectional prevalence and recency of infection tests as described by Kassanjee, et al. "A new general biomarker-based incidence estimator," \emph{Epidemiology} (2012). Function parameters must be specified to include assay test characteristics and survey results as proportions. Confidence intervals are computed via Delta method approximation, except when Boot=TRUE is specified, in which case confidence intervals are generated by empirical bootstrap resampling. Inputs must be in appropriate ranges for appropriate units. Extreme input values may make calculation impossible, and if entered will elicit error notices.
#' @return Incidence estimate, confidence interval, relative standard error of estimate and of assay characteristics MDRI and FRR. If multiple surveys are entered, function returns said results, as well as estimates of incidence  differences, confidence intervals of differences, difference relative standard errors, and p-values testing the hypothesis that the difference in incidence measures are zero. Theoretical relative standard error of incidence and incidence difference at infinite sample size is returned only if Boot=FALSE, as that calculation relies on the asymptotic behavior of components of the Delta method approximation, and is not calculable from bootstrapped values.
#' @examples
#' incBYcounts(N = c(5000) ,N_H = 1000, N_testR = 1000, N_R = 70,
#' BSoot = FALSE, BMest = "MDRI.FRR.indep", MDRI = 200, RSE_MDRI = 0.05,
#' FRR = 0.01, RSE_FRR = 0.2, BigT = 730)
#'
#'
#' incBYcounts(N = c(4000,4000,4050) ,N_H = c(1010,1000,900),
#' N_testR = c(1000,1000,880), N_R = c(60,70,50), Boot=TRUE,
#' BMest="same.test", MDRI = 210, RSE_MDRI = 0.05, FRR = 0.005,
#' RSE_FRR = 0.19, BigT = 700)
#'
#'
#' incBYcounts(N = c(4000,4000,4000) ,N_H = c(1050,1090),
#' N_testR = c(1000,1000), N_R = c(60,67), Boot=FALSE, BMest="FRR.indep",
#' MDRI = 220, RSE_MDRI = 0.05, FRR = c(0.005,0.005,0.01), RSE_FRR = 0.19,
#' BigT = 610)
#'
#' @export
incBYcounts<-function(N, N_H, N_testR, N_R,
                      DE_H = 1, DE_R = 1,
                      BS_Count = 10000, Boot = FALSE, alpha = 0.05,
                      BMest = "same.test", MDRI, RSE_MDRI, FRR, RSE_FRR,
                      BigT = 730, Covar_HR = 0){

  if(sum(BMest==c("same.test", "FRR.indep", "MDRI.FRR.idep"))==0){
    stop("BMest option must be same.test, FRR.indep, or MDRI.FRR.idep")
  }

  counts.to.prev<-prevBYcounts(N=N, N_H=N_H, N_testR=N_testR,N_R=N_R, DE_H=DE_H, DE_R=DE_R)

  recencyI(BS_Count=BS_Count, Boot=Boot, alpha=alpha, BMest=BMest,
           PrevH=counts.to.prev$PrevH, RSE_PrevH=counts.to.prev$RSE_PrevH,
           PrevR=counts.to.prev$PrevR, RSE_PrevR=counts.to.prev$RSE_PrevR,
           MDRI=MDRI, RSE_MDRI=RSE_MDRI, FRR=FRR, RSE_FRR=RSE_FRR,
           BigT=BigT, Covar_HR=Covar_HR)
}







# Begin Examples ----------------------------------------------------------

######## -- The rest of this code gives example to show the function agrees with spreadsheet answers,
#and to see when and how the function breaks, what error messages it outputs.
#Final R package must omit this code


probs<-prevBYcounts (N=c(5000,5000), N_H=c(1000,1000), N_testR=c(1000,900), N_R=c(100,70))

################### == Call - DM only ==######################################################
recencyI  (BS_Count=10000,
           Boot=FALSE,
           BMest="same.test",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)

#(single survey)
recencyI  (BS_Count=10000,
           Boot=FALSE,
           BMest="same.test",
           PrevH=probs[1,1], RSE_PrevH=probs[1,3],
           PrevR=probs[1,2], RSE_PrevR=probs[1,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)



##############################################################################################

################### == Call - BS only ==######################################################
recencyI  (BS_Count=10000,
           Boot=TRUE,
           BMest="same.test",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)

#boot with a larger sample
recencyI  (BS_Count=100000,
           Boot=TRUE,
           BMest="same.test",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730, Covar_HR = 0.00002)

#boot with a larger sample and different testing scheme
recencyI  (BS_Count=100000,
           Boot=TRUE,
           BMest="FRR.indep",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=c(0.01,0.008), RSE_FRR=c(0.2,0.21),
           BigT=730, Covar_HR = 0.00002)
##############################################################################################



########== Check with spread sheets """"""same.test, three surveys""""""==###########
probs<-prevBYcounts(N=c(5000,5000,3000), N_H=c(1000,1000,1000), N_testR=c(1000,1000,900), N_R=c(100,70,120))

recencyI  (BS_Count=10000,
           Boot=FALSE,
           BMest="same.test",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)

########== Check with spread sheets """"""FRR.indep""""""==###########
probs<-prevBYcounts (N=c(5000,5000,3000), N_H=c(1000,1000,1000), N_testR=c(1000,1000,900), N_R=c(100,70,120))

recencyI  (BS_Count=10000,
           Boot=FALSE,
           BMest="FRR.indep",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)

########== Check with spread sheets """"""MDRI.FRR.indep""""""==###########
probs<-prevBYcounts (N=c(5000,5000,3000), N_H=c(1000,1000,1000), N_testR=c(1000,1000,900), N_R=c(100,70,120))

recencyI  (BS_Count=10000,
           Boot=FALSE,
           BMest="MDRI.FRR.indep",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=c(0.01,0.02,0.001), RSE_FRR=0.2,
           BigT=730)


#######== bootstrap with warining from FRR=0 ==###########
probs<-prevBYcounts (N=c(5000,5000,3000), N_H=c(1000,1000,1000), N_testR=c(1000,1000,900), N_R=c(100,70,120))

recencyI  (BS_Count=10000,
           Boot=TRUE,
           BMest="MDRI.FRR.indep",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0.01, RSE_FRR=0.2,
           BigT=730)

recencyI  (BS_Count=10000,
           Boot=FALSE,
           BMest="same.test",
           PrevH=probs[,1], RSE_PrevH=probs[,3],
           PrevR=probs[,2], RSE_PrevR=probs[,4],
           MDRI=200, RSE_MDRI=0.05,
           FRR=0, RSE_FRR=0,
           BigT=730)


########== incidence by counts, single survey ==###########
incBYcounts (N=5000, N_H=1000, N_testR=1000, N_R=100,
             DE_H=1, DE_R=1,
             BS_Count=10000, Boot=FALSE,
             BMest="same.test", MDRI=200, RSE_MDRI=.05, FRR=0.01, RSE_FRR=0.2,
             BigT=730, Covar_HR=0)

########== incidence by counts, two surveys ==###########
incBYcounts (N=c(5000,5000), N_H=c(1000,1000), N_testR=c(1000,1000), N_R=c(100,70),
             DE_H=1, DE_R=1,
             BS_Count=10000, Boot= FALSE,
             BMest="same.test", MDRI=200, RSE_MDRI=.05, FRR=0.01, RSE_FRR=0.06,
             BigT=730, Covar_HR=0)

########== incidence by counts, two surveys, FRR independent ==###########
incBYcounts (N=c(5000,5000), N_H=c(1000,1000), N_testR=c(1000,950), N_R=c(100,70),
             DE_H=1, DE_R=1,
             BS_Count=10000, Boot= FALSE,
             BMest="FRR.indep", MDRI=200, RSE_MDRI=.05, FRR=0, RSE_FRR=0.05,
             BigT=730, Covar_HR=0)



# End Examples ------------------------------------------------------------










