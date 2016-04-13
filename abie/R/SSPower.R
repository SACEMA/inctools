#SCRIPT UNFINISHED

#This script is a function that takes either a sample size and gives power for a test of two incidences,
#or takes power and gives the required sample size necessary.





#add N_testR to calculation...



#' Power and Sample Size calculation for assay-based incidence estimation
#'
#' @param I1 Predicted incidence of HIV in survey 1
#' @param I2 Predicted incidence of HIV in survey 2
#' @param PrevH1 Predicted prevalence of HIV in survey 1
#' @param PrevH2 Predicted prevalence of HIV in survey 2
#' @param n1 Sample size for survey 1
#' @param n2 Sample size for survey 2
#' @param alpha Significance level for test (default alpha=0.05)
#' @param Power Desired power used to caluculate a sample size for the surveys. Default is 0.80, meaning the function outputs the necessary sample size to achieve stated power for a test of differences in incidence.
#' @param SS Sample size. Default is NULL, meaning the function takes a power argument and outputs sample size needed to achieve power level for test of differences for incidence. If parameter is set to a number, then function will return the power achieved by the test of differences under that sample size.
#' @param DE_H Design effect of HIV-prevalence test (vector/integer)
#' @param DE_R Design effect of recency test (vector/integer)
#' @param BMest Biomarker estimation by one the 3 options "same.test"(=default), "FRR.indep", "MDRI.FRR.indep" (string)
#' @param MDRI mean duration of recent infection [days] (vector/integer)
#' @param RSE_MDRI Relative standard error of MDRI [days] (vector/integer)
#' @param FRR False recent rate (vector/integer)
#' @param RSE_FRR Relative standard error of FRR (vector/integer)
#' @param BigT post-infection time cut-off true vs false recent [days] default 730 days (integer)
#' @return Sample size necessary for a given power level for testing the hypothesis that the incidence rates are identical between populations; alternatively, the power of said test under a particular sample size scenario. Function also returns implied statistics from input values, confidence limits, and population counts.
#' @details ...
#'
#' @examples
#'SSPower(I1=0.05, I2=0.03, PrevH1=0.20, PrevH2=0.20, n1=5000,
#'n2=5000, alpha=0.05, Power="out", SS=NULL, DE_H=1, DE_R=1,
#'BMest="same.test", MDRI=200, RSE_MDRI=0.05, FRR=0.01,
#'RSE_FRR=0.20, BigT=730)
#'
#'
#'SSPower(I1=0.05, I2=0.03, PrevH1=0.20, PrevH2=0.20, alpha=0.05,
#'Power=0.80, SS="out", DE_H=1, DE_R=1, BMest="FRR.indep",
#'MDRI=200, RSE_MDRI=0.05, FRR=0.01, RSE_FRR=0.20, BigT=730)
#' @export


SSPower <- function (I1, I2, PrevH1, PrevH2, n1, n2, alpha=0.05, Power=0.80, SS="out", DE_H=1, DE_R=1,
                      BMest="same.test", MDRI, RSE_MDRI, FRR, RSE_FRR,
                      BigT=730){

  # stopifnot (PrevH<=1     & PrevH>=0)
  # stopifnot (PrevR<=1     & PrevR>=0)
  # stopifnot (RSE_PrevH<=1 & RSE_PrevH>=0)
  # stopifnot (RSE_PrevR<=1 & RSE_PrevR>=0)
  # stopifnot (MDRI>=0)
  # stopifnot (RSE_MDRI<=1  & RSE_MDRI>=0)
  # stopifnot (FRR<=1       & FRR>=0)
  # stopifnot (RSE_FRR<=1   & RSE_FRR>=0)

  if(sum(BMest==c("same.test", "FRR.indep", "MDRI.FRR.indep"))==0){
    stop("BMest option must be same.test, FRR.indep, or MDRI.FRR.idep")
  }


  if (length(MDRI)>length(FRR)) {stop("number of inputs for MDRI is larger than number of inputs for FRR")}

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

  if (BigT<=182)                {warning ("BigT is smaller than half a year")}

  no_s <- 2 #dimension of inputs (number of surveys)
  if (length(MDRI)==1)     {MDRI <- rep(MDRI, times=no_s)}         else {MDRI=MDRI}
  if (length(FRR)==1)      {FRR  <- rep(FRR, times=no_s)}          else {FRR=FRR}
  if (length(RSE_MDRI)==1) {RSE_MDRI <- rep(RSE_MDRI, times=no_s)} else {RSE_MDRI=RSE_MDRI}
  if (length(RSE_FRR)==1)  {RSE_FRR  <- rep(RSE_FRR, times=no_s)}  else {RSE_FRR=RSE_FRR}
  if (length(DE_H)==1)  {DE_H  <- rep(DE_H, times=no_s)}  else {DE_H=DE_H}
  if (length(DE_R)==1)  {DE_R  <- rep(DE_R, times=no_s)}  else {DE_R=DE_R}
  if (length(CR)==1)  {CR  <- rep(CR, times=no_s)}  else {CR=CR}


  stopifnot(no_s==length(MDRI)    & no_s==length(RSE_MDRI)  & no_s==length(FRR)   &
              no_s==length(RSE_FRR)  & length(BigT)==1)


  MDRI<-MDRI/365.25
  BigT<-BigT/365.25
  if(Power=="out")N <- c(n1,n2)
  PrevH <- c(PrevH1,PrevH2)
  I<-c(I1,I2)
  deltaI_Est <- I[1]-I[2]

  HIV.neg<- 1-PrevH
  PrevR <- I*(1-PrevH)*(MDRI-FRR*BigT)+(FRR*PrevH)

  DM_Var_MDRI  <- (MDRI*RSE_MDRI)^2
  DM_Var_FRR   <- (FRR*RSE_FRR)^2

  MDRI.CI <- 365.25*data.frame(Implied.MDRI.CI.low=qnorm(alpha/2, mean=MDRI, sd=sqrt(DM_Var_MDRI)),Implied.MDRI.CI.up=qnorm(1-alpha/2, mean=MDRI, sd=sqrt(DM_Var_MDRI)))
  FRR.CI <- data.frame(Implied.FRR.CI.low=qnorm(alpha/2, mean=FRR, sd=sqrt(DM_Var_FRR)),Implied.FRR.CI.up=qnorm(1-alpha/2, mean=FRR, sd=sqrt(DM_Var_FRR)))

#
#   DM_Var_PrevH<-PrevH*(1-PrevH)*DE_H/N
#   DM_Var_PrevR<-PrevR*(1-PrevR)*DE_R/N


# #
#     #next few lines make delta method matrix
#     fot_Mat  <- matrix (nrow=no_s, ncol=4)
#     DM_Var_I <- vector(length=no_s)
#     for (i in c(1:no_s)) {
#       fot_Mat[i,] <- DM_FirstOrderTerms (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i], frr=FRR[i], bigt=BigT)
#       #DM_FirstOrderTerms gives fot_prevH, fot_prevR, fot_mdri, fot_frr, so FOT for each prev. survey input
#       #and each row of fot_Mat has the FOT for prev. survey i.
#       DM_Var_I[i] <- (fot_Mat[i,1]^2)*DM_Var_PrevH[i] + (fot_Mat[i,2]^2)*DM_Var_PrevR[i] +
#         (fot_Mat[i,3]^2)*DM_Var_MDRI[i]  + (fot_Mat[i,4]^2)*DM_Var_FRR[i]
#     }
#     #I tested this function explicitely. It works, and works in the other function scripts as well.
#

# #     DM_Var_deltaI <- DM_VAR_deltaI (BMest=BMest, fot_prevH1=fot_Mat[1,1], fot_prevH2=fot_Mat[2,1],
# #                                                      fot_prevR1=fot_Mat[1,2],   fot_prevR2=fot_Mat[2,2],
# #                                                      fot_mdri1=fot_Mat[1,3],    fot_mdri2=fot_Mat[2,3],
# #                                                      fot_frr1=fot_Mat[1,4],     fot_frr2=fot_Mat[2,4],
# #                                                      dm_var_prevH1=DM_Var_PrevH[1], dm_var_prevH2=DM_Var_PrevH[2],
# #                                                      dm_var_prevR1=DM_Var_PrevR[1], dm_var_prevR2=DM_Var_PrevR[2],
# #                                                      dm_var_mdri1=DM_Var_MDRI[1],   dm_var_mdri2=DM_Var_MDRI[2],
# #                                                      dm_var_frr1=DM_Var_FRR[1],     dm_var_frr2=DM_Var_FRR[2])
# #
# #
# #     DM_Var_deltaI.infSS <- DM_VAR_deltaI.infSS (BMest=BMest,
# #                                     fot_mdri1=fot_Mat[1,3],    fot_mdri2=fot_Mat[2,3],
# #                                     fot_frr1=fot_Mat[1,4],     fot_frr2=fot_Mat[2,4],
# #                                     dm_var_mdri1=DM_Var_MDRI[1],   dm_var_mdri2=DM_Var_MDRI[2],
# #                                     dm_var_frr1=DM_Var_FRR[1],     dm_var_frr2=DM_Var_FRR[2])
#
# #
# #   DM_SD_I      <- sqrt(DM_Var_I)
# #   DM_SD_deltaI <- sqrt(DM_Var_deltaI)
# #   DM_SD_deltaI.infSS <- sqrt(DM_Var_deltaI.infSS)
# #   Var_I        <-  DM_Var_I
# #   RSE_I        <- sqrt(Var_I)/I
# #   SD_I         <- sqrt(Var_I)
# #   Var_deltaI   <- DM_Var_deltaI #+  BS_Var_deltaI
# #   RSE_deltaI   <- sqrt(Var_deltaI)/abs(deltaI_Est) #sqrt(Var_deltaI)/abs(deltaI_Est_Vec)
# #   RSE_deltaI.infSS   <- sqrt(DM_Var_deltaI.infSS)/abs(deltaI_Est)
# #   SD_deltaI    <- sqrt(Var_deltaI)
#



#Instead of trying to get the fot function to work, I explicitely wrote out the commands for the Var[I], Var[D-I]
if(Power=="out"){
  Var_I <- I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))
                + (RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2
                +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2
                )
  RSE_I <- sqrt(Var_I)/I
  Var_delta_I <-
    if(BMest=="MDRI.FRR.indep"){Var_I[1]+Var_I[2]
    } else if(BMest=="FRR.indep") { sum(I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))))  +   ((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2  +    sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)
        } else if(BMest=="same.test"){ sum(I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))))  +   ((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)   +   ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2) }

  RSE_deltaI <- sqrt(Var_delta_I)/abs(deltaI_Est)
  RSE_deltaI.infSS <-  if(BMest=="MDRI.FRR.indep"){sqrt(((RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2 +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)*I^2)/I
  } else if(BMest=="FRR.indep") {   sqrt(((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2  + sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2))/deltaI_Est
  } else if(BMest=="same.test"){    ((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)   +   ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2) }



  CI.low <- qnorm(alpha/2, mean=I, sd=sqrt(Var_I))
  CI.up <- qnorm(1-alpha/2, mean=I, sd=sqrt(Var_I))

  deltaI_CI<-NULL
  deltaI_CI[1] <- qnorm(alpha/2, mean=deltaI_Est, sd=sqrt(Var_delta_I))
  deltaI_CI[2] <- qnorm(1-alpha/2, mean=deltaI_Est, sd=sqrt(Var_delta_I))

  ss.power <- 1-pnorm(q=qnorm(1-alpha/2), mean=1/RSE_deltaI, sd=1)
  Power.infSS <-1-pnorm(q=qnorm(1-alpha/2), mean=1/RSE_deltaI.infSS, sd=1)

  #if(ss.power<0.7){warning("Probability of correct inference less than 70%")}
  #this warning is in spreadsheets, but I don't like it. User expected to have an idea about power...


  if(BMest=="FRR.indep"){
  output <- list(Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3))),
                 Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), Implied.RSE_I=round(RSE_I,3), Implied.CI.low=round(CI.low,3), Implied.CI.up=round(CI.up,3)),
                 Implied.FRR.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), Implied.MDRI.CI.low=round(MDRI.CI[1,1],3), Implied.MDRI.CI.up=round(MDRI.CI[1,2],3)),
                 Implied.MDRI.Statistics=data.frame(Given.FRR=round(FRR,3), Implied.FRR.CI.low=round(FRR.CI[,1],3), Implied.FRR.CI.up=round(FRR.CI[,2],3)),
                 Implied.Subject.Counts=data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) )
} else if(BMest=="same.test") {output <- list(Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3))),
                      Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), Implied.RSE_I=round(RSE_I,3), Implied.CI.low=round(CI.low,3), Implied.CI.up=round(CI.up,3)),
                      Implied.FRR.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), Implied.MDRI.CI.low=round(MDRI.CI[1,1],3), Implied.MDRI.CI.up=round(MDRI.CI[1,2],3)),
                      Implied.MDRI.Statistics=data.frame(Given.FRR=round(FRR[1],3), Implied.FRR.CI.low=round(FRR.CI[1,1],3), Implied.FRR.CI.up=round(FRR.CI[1,2],3)),
                      Implied.Subject.Counts=data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) )
} else if(BMest=="MDRI.FRR.indep"){ output <- list(Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3))),
                                                  Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), Implied.RSE_I=round(RSE_I,3), Implied.CI.low=round(CI.low,3), Implied.CI.up=round(CI.up,3)),
                                                  Implied.FRR.Statistics=data.frame(Given.MDRI=round(MDRI*365.25,3), Implied.MDRI.CI.low=round(MDRI.CI[,1],3), Implied.MDRI.CI.up=round(MDRI.CI[,2],3)),
                                                  Implied.MDRI.Statistics=data.frame(Given.FRR=round(FRR,3), Implied.FRR.CI.low=round(FRR.CI[,1],3), Implied.FRR.CI.up=round(FRR.CI[,2],3)),
                                                  Implied.Subject.Counts=data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) )
}
} else if(SS=="out"){
if (BMest=="same.test")
  {SS <- ceiling(sum(I^2*((1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))))/
                    ((.02/(qnorm(1-0.8)-qnorm(1-alpha/2)) )^2 -
                    ((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)  -
                    ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2)
                     ))
} else if(BMest=="FRR.indep")
  { SS <- ceiling(sum(I^2*((1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2)))) /
                  ((deltaI_Est/(qnorm(1-0.858)-qnorm(1-alpha/2)))^2  -  ((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2 -
                  sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)
                  ))
} else if(BMest=="MDRI.FRR.indep")
  { SS <- ceiling(sum(I^2*((1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2)))) /
                      ((deltaI_Est/(qnorm(1-0.858)-qnorm(1-alpha/2)))^2  -  sum(I^2*(RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2) -
                      sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2 )
                      ))
}
  #Now based on derived common SS, output implied summary statistics
  N<-c(SS,SS)
  Var_I <- I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))
                + (RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2
                +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2
  )
  RSE_I <- sqrt(Var_I)/I
  Var_delta_I <-
    if(BMest=="MDRI.FRR.indep"){Var_I[1]+Var_I[2]
    } else if(BMest=="FRR.indep") { sum(I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))))  +   ((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2  +    sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)
    } else if(BMest=="same.test"){ sum(I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))))  +   ((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)   +   ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2) }

  RSE_deltaI <- sqrt(Var_delta_I)/abs(deltaI_Est)
  RSE_deltaI.infSS <-  if(BMest=="MDRI.FRR.indep"){ sqrt(((RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2 +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)*I^2)/I
  } else if(BMest=="FRR.indep") {   sqrt(((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2  + sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2))/deltaI_Est
  } else if(BMest=="same.test"){    sqrt(((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)   +   ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2))/deltaI_Est }

  CI.low <- qnorm(alpha/2, mean=I, sd=sqrt(Var_I))
  CI.up <- qnorm(1-alpha/2, mean=I, sd=sqrt(Var_I))

  deltaI_CI<-NULL
  deltaI_CI[1] <- qnorm(alpha/2, mean=deltaI_Est, sd=sqrt(Var_delta_I))
  deltaI_CI[2] <- qnorm(1-alpha/2, mean=deltaI_Est, sd=sqrt(Var_delta_I))

  ss.power <- 1-pnorm(q=qnorm(1-alpha/2), mean=1/RSE_deltaI, sd=1)
  Power.infSS <-1-pnorm(q=qnorm(1-alpha/2), mean=1/RSE_deltaI.infSS, sd=1)


  if(BMest=="FRR.indep"){
    output <- list(Minimum.Common.SS=SS,
                   Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3))),
                   Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), Implied.RSE_I=round(RSE_I,3), Implied.CI.low=round(CI.low,3), Implied.CI.up=round(CI.up,3)),
                   Implied.FRR.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), Implied.MDRI.CI.low=round(MDRI.CI[1,1],3), Implied.MDRI.CI.up=round(MDRI.CI[1,2],3)),
                   Implied.MDRI.Statistics=data.frame(Given.FRR=round(FRR,3), Implied.FRR.CI.low=round(FRR.CI[,1],3), Implied.FRR.CI.up=round(FRR.CI[,2],3)),
                   Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
  } else if(BMest=="same.test") {output <- list(Minimum.Common.SS=SS,
                                                Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3))),
                                                Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), Implied.RSE_I=round(RSE_I,3), Implied.CI.low=round(CI.low,3), Implied.CI.up=round(CI.up,3)),
                                                Implied.FRR.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), Implied.MDRI.CI.low=round(MDRI.CI[1,1],3), Implied.MDRI.CI.up=round(MDRI.CI[1,2],3)),
                                                Implied.MDRI.Statistics=data.frame(Given.FRR=round(FRR[1],3), Implied.FRR.CI.low=round(FRR.CI[1,1],3), Implied.FRR.CI.up=round(FRR.CI[1,2],3)),
                                                Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
  } else if(BMest=="MDRI.FRR.indep"){ output <- list(Minimum.Common.SS=SS,
                                                     Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3))),
                                                     Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), Implied.RSE_I=round(RSE_I,3), Implied.CI.low=round(CI.low,3), Implied.CI.up=round(CI.up,3)),
                                                     Implied.FRR.Statistics=data.frame(Given.MDRI=round(MDRI*365.25,3), Implied.MDRI.CI.low=round(MDRI.CI[,1],3), Implied.MDRI.CI.up=round(MDRI.CI[,2],3)),
                                                     Implied.MDRI.Statistics=data.frame(Given.FRR=round(FRR,3), Implied.FRR.CI.low=round(FRR.CI[,1],3), Implied.FRR.CI.up=round(FRR.CI[,2],3)),
                                                     Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
  }

   }




  return(output)
}








##################### ---  Test values against spreadsheets (Find Power)---- #######################
I1 <- 0.05
I2 <- 0.03
PrevH1 <- .2
PrevH2 <- .2
Power="out"
MDRI     <- 200
RSE_MDRI <- 0.05
FRR      <- 0.01
RSE_FRR  <- 0.2
BigT     <- 730
CR       <- 1
DE_H     <- 1
DE_R     <- 1
n1       <- 5000
n2       <-5000
alpha=0.05
BMest="same.test"

SSPower(I1, I2, PrevH1, PrevH2, n1, n2, alpha=0.05, Power="out", SS=NULL, DE_H=1, DE_R=1,
                     BMest="same.test", MDRI, RSE_MDRI, FRR, RSE_FRR,
                     BigT=730)


my.data<-SSPower(I1, I2, PrevH1, PrevH2, n1, n2, alpha=0.05, Power="out", SS=NULL, DE_H=1, DE_R=1,
                 BMest="same.test", MDRI, RSE_MDRI, FRR, RSE_FRR,
                 BigT=730)
my.data$Inc.Difference.Statistics
my.data[[1]]
my.data[[1]][4]
my.data[[1]]$Power
names(my.data)
str(my.data)



##################### ---  Test values against spreadsheets (Find Power)---- #######################

I1 <- 0.05
I2 <- 0.03
PrevH1 <- .2
PrevH2 <- .15
Power= .8
SS = "out"
MDRI     <- 200
RSE_MDRI <- 0.05
FRR      <- 0.01
RSE_FRR  <- 0.2
BigT     <- 730
CR       <- 1
DE_H     <- 1
DE_R     <- 1
n1       <- NULL
n2       <-NULL
alpha=0.05
BMest="same.test"



test.SS <- SSPower(I1, I2, PrevH1, PrevH2, n1, n2, alpha=0.05, Power=.8, SS="out", DE_H=DE_H, DE_R=DE_R,
                 BMest="same.test", MDRI, RSE_MDRI, FRR, RSE_FRR,
                 BigT=730)

test.SS






##################### ---  Test values against spreadsheets (Find Sample Size)---- #######################












