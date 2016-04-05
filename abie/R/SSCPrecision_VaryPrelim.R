####== Sample Size calculation for desired precision of I estimation
#######################################################################################################
####== PrevH:    prevalence of HIV
####== PrevR:    prevalence of recency when HIV positive and therefore tested for recency
####== I:        incidence expected
####== CR:       coverage rate; probability (0 - 1) of being tested for recency when positive
####==
####== Covar_HR: covariance of probs being postive & being recent of survey 1&2 (vector)
####== MDRI:     mean duration of recent infection [days] (vector/integer)
####== RSE_MDRI: Relative standard error of MDRI   [days] (vector/integer)
####== FRR:      False recent rate (vector/integer)
####== RSE_FRR:  Relative standard error of FRR (vector/integer)
####== BigT:     post-infection time cut-off true vs false recent [days] default 730 days (integer)
####== DE_H:     Design effect of HIV-prevalence test (vector/integer)
####== DE_R:     Design effect of recency test (vector/integer)
####== BS_Vars:  Variables to be bootstrapped e.g.c("PrevH", "PrevR", "MDRI", "FRR") (string vector)
####==
#######################################################################################################

#This function is to be separate from SS-Power function. This function takes as argument test characteristics,
#and returns EITHER sample size for a given precision, or a given precision for a given sample size.

MDRI     <- c(200,250)
RSE_MDRI <- 0.05
FRR      <- 0.01
RSE_FRR  <- 0.2
BigT     <- 730
I        <- 0.015
RSE_I    <- 0.25
PrevH    <- c(0.1,0.5)
CR       <- 1
DE_H     <- 1
DE_R     <- 1
n        <- "out"
step <- 5

#Note to self: 'step' is number of steps between minimum I and maximum I in the calculation of
#a range of output. So supply a vector or max/min theoretical incidences, and the function gives
#a range of values (step number of values) for the output. This can be done for all input variables.
#DO WE WANT THIS OPTION??

#Only gives option of n or RSE_I to be 'out' or not. So only gives precision or n. Makes sense.
#Either give function n to get RSE_I or give it RSE_I to get n.

#This function covers what was made in spreadsheets ABIE_v3_Test_Performance_Calculator (which takes SS
# and other variables of test, hypothetical data, and gives precision), and




#' Sample Size or Precision Calculation
#'
#' @param I Expected Incidence
#' @param RSE_I Relative Standard Error of Incidence Estimate
#' @param PrevH Prevalence of HIV
#' @param CR Coverage rate: probability (0-1) of being tested for recency when positive for HIV
#' @param MDRI mean duration of recent infection [days] (vector/integer)
#' @param RSE_MDRI Relative standard error of MDRI   [days] (vector/integer)
#' @param FRR False recent rate (vector/integer)
#' @param RSE_FRR Relative standard error of FRR (vector/integer)
#' @param BigT post-infection time cut-off true vs false recent [days] default 730 days (integer)
#' @param DE_H Design effect of HIV-prevalence test (vector/integer)
#' @param DE_R Design effect of recency test (vector/integer)
#' @param n Sample Size, either given hypothetical value or to be determined by function
#' @param step number of steps between minimum I and maximum I in the calculation of a range of output
#' @return Either sample size necessary for a given precision under a given set of testing characteristics and a hypothetical prevalence/incidence scenario, or precision under a particular sample size scenario, with a given hypothetical prevalence/incidence scenario.
#' @details
#'
#' Summarizes performance of a recent infection test (into a standard error of the incidence estimate), given estimated test properties and the prevalence/incidence in a hypothetical context; or gives sample size necessary for a given level of estimator precision
#'
#' @examples
#' SSCprecision(I = 0.015, RSE_I = 0.25, PrevH = 0.2, CR = 1, MDRI = 200, RSE_MDRI = 0.05,
#' FRR = 0.01, RSE_FRR = 0.2, BigT = 730, DE_H = 1, DE_R = 1, n = "out", step = 5)
#' @export
#'
SSCprecision <- function ( I              ,
                           RSE_I          ,
                           PrevH          ,
                           CR             ,
                           MDRI           ,
                           RSE_MDRI       ,
                           FRR            ,
                           RSE_FRR        ,
                           BigT           = 730,
                           DE_H           = 1,
                           DE_R           = 1,
                           n              = "out",
                           step           = 5)
{
  test <- c( I,RSE_I,PrevH,CR ,MDRI ,RSE_MDRI,FRR,RSE_FRR,BigT,DE_H,DE_R,n,step)
  if (length(test) > 15) {
    stop("only a maximum of 2 variables are allowed to vary")
  }

  if (length(I)>2        | length(I)<1)       {("specifiy (only) min & max values for I")}
  if (length(RSE_I)>2    | length(RSE_I)<1)   {("specifiy (only) min & max values for RSE_I")}
  if (length(PrevH)>2    | length(PrevH)<1)   {("specifiy (only) min & max values for PrevH")}
  if (length(CR)>2       | length(CR)<1)      {("specifiy (only) min & max values for CR")}
  if (length(MDRI)>2     | length(MDRI)<1)    {("specifiy (only) min & max values for MDRI")}
  if (length(RSE_MDRI)>2 | length(RSE_MDRI)<1){("specifiy (only) min & max values for RSE_MDRI")}
  if (length(FRR)>2      | length(FRR)<1)     {("specifiy (only) min & max values for FRR")}
  if (length(RSE_FRR)>2  | length(RSE_FRR)<1) {("specifiy (only) min & max values for RSE_FRR")}
  if (length(BigT)>2     | length(BigT)<1)    {("specifiy (only) min & max values for BigT")}
  if (length(DE_H)>2     | length(DE_H)<1)    {("specifiy (only) min & max values for DE_H")}
  if (length(DE_R)>2     | length(DE_R)<1)    {("specifiy (only) min & max values for DE_R")}
  if (length(n)>2        | length(n)<1)       {("specifiy (only) min & max values for n")}

#THIS SECTION STARTING HERE NEEDS UPDATED ERROR NOTES
  if (is.numeric(I)>0)        {stopifnot (I<=1        & I>=0)}
  if (is.numeric(RSE_I)>0)    {stopifnot (RSE_I<=1    & RSE_I>=0)}
  if (is.numeric(PrevH)>0)    {stopifnot (PrevH<=1    & PrevH>=0)}
  if (is.numeric(CR)>0)       {stopifnot (CR<=1       & CR>=0)}
  if (is.numeric(MDRI)>0)     {stopifnot (MDRI>=0)}
  if (is.numeric(RSE_MDRI)>0) {stopifnot (RSE_MDRI<=1 & RSE_MDRI>=0)}
  if (is.numeric(FRR)>0)      {stopifnot (FRR<=1      & FRR>=0)}
  if (is.numeric(RSE_FRR)>0)  {stopifnot (RSE_FRR<=1  & RSE_FRR>=0)}
  if (is.numeric(DE_H)>0)     {stopifnot (DE_H>=1)}
  if (is.numeric(DE_R)>0)     {stopifnot (DE_R>=1)}
  if (is.numeric(n)>0)        {stopifnot (n>100)}

#CAN BE MADE MORE SUCCINCT
    if (length(BigT)==1) {
    if (BigT<=182)          {
      warning ("BigT is smaller than half a year")
      }
  } else {
    for (i in c(1,2)) {
      if (BigT[i]<=182)          {
        warning ("BigT is smaller than half a year")
      }
    }
  }

  vary1 <- NULL
  vary2 <- NULL
  if (length(I)==2)      {
    I <- matrix(rep(seq(from=min(I),to=max(I),length.out=step),times=step),ncol=step,nrow=step)
    if (length(vary1)==0)         {
      vary1 <- I
      vary_name1 <- "I"
    } else {
      vary2 = I = t(I)
      vary_name2 = "I"
    }
  }
  if (length(RSE_I)==2)   {
    RSE_I <- matrix(rep(seq(from=min(RSE_I),to=max(RSE_I),length.out=step),times=step),ncol=step,nrow=step)
    if (length(vary1)==0){
      vary1 <- RSE_I
      vary_name1 <- "I"
    } else {
      vary2 = RSE_I = t(RSE_I)
      vary_name2 = "I"
    }
  }
  if (length(PrevH)==2)   {
    PrevH  <- matrix(rep(seq(from=min(PrevH),to=max(PrevH),length.out=step),times=step),ncol=step,nrow=step)
    if (length(vary1)==0){
      vary1 <- PrevH
      vary_name1 <- "PrevH"
    } else {
      vary2 = PervH = t(PrevH)
      vary_name2 <- "PrevH"
    }
  }
  if (length(CR)==2) {
    CR <- matrix(rep(seq(from=min(CR),to=max(CR),length.out=step),times=step),ncol=step,nrow=step)
    if (length(vary1)==0)         {
      vary1 <- CR
      vary_name1 <- "CR"
    } else {
      vary2 = CR = t(CR)
      vary_name2 = "CR"
    }
  }
  if (length(MDRI)==2) {
    MDRI <- matrix(rep(seq(from=min(MDRI),to=max(MDRI),length.out=step),times=step),ncol=step,nrow=step)
    if (length(vary1)==0)         {
      vary1 <- MDRI
      vary_name1 <- "MDRI"
    } else {
      vary2 = MDRI = t(MDRI)
      vary_name2 = "MDRI"
    }
  }
  if (length(RSE_MDRI)==2){
    RSE_MDRI <- matrix(rep(seq(from=min(RSE_MDRI),to=max(RSE_MDRI),length.out=step),times=step),ncol=step,nrow=step)
    if (length(vary1)==0)         {
      vary1 <- RSE_MDRI
      vary_name1 <- "RSE_MDRI"
    } else {
      vary2 = RSE_MDRI = t(RSE_MDRI)
      vary_name2 = "RSE_MDRI"
    }
  }
  if (length(FRR)==2)  {
    FRR <- matrix(rep(seq(from=min(FRR),to=max(FRR),length.out=step),times=step),nrow=step,ncol=step)
    if (length(vary1)==0)         {
      vary1 <- FRR
      vary_name1 <- "FRR"
    } else {
      vary2 = FRR = t(FRR)
      vary_name2 = "FRR"
    }
  }
  if (length(RSE_FRR)==2) {
    RSE_FRR <- matrix(rep(seq(from=min(RSE_FRR),to=max(RSE_FRR),length.out=step),times=step),nrow=step,ncol=step)
    if (length(vary1)==0)         {
      vary1 <- RSE_FRR
      vary_name1 <- "RSE_FRR"
    } else {
      vary2 = RSE_FRR = t(RSE_FRR)
      vary_name2 = "RSE_FRR"
    }
  }
  if (length(BigT)==2)  {
    BigT <- matrix(rep(seq(from=min(BigT),to=max(BigT),length.out=step),times=step),nrow=step,ncol=step)
    if (length(vary1)==0)         {
      vary1 <- BigT
      vary_name1 <- "BigT"
    } else {
      vary2 = BigT = t(BigT)
      vary_name2 = "BigT"
    }
  }
  if (length(DE_H)==2)    {
    DE_H <- matrix(rep(seq(from=min(DE_H),to=max(DE_H),length.out=step),times=step),nrow=step,ncol=step)
    if (length(vary1)==0)         {
      vary1 <- DE_H
      vary_name1 <- "DE_H"
    } else {
      vary2 = DE_H = t(DE_H)
      vary_name2 = "DE_H"
    }
  }
  if (length(DE_R)==2)  {
    DE_R <- matrix(rep(seq(from=min(DE_R),to=max(DE_R),length.out=step),times=step),nrow=step,ncol=step)
    if (length(vary1)==0)        {
      vary1 <- DE_R
      vary_name1 <- "DE_R"
    } else {
      vary2 = DE_R = t(DE_R)
      vary_name2 = "DE_R"
    }
  }
  if (length(n)==2)   {
    n <- matrix(rep(seq(from=min(n),to=max(n),length.out=step),times=step),nrow=step,ncol=step)
    if (length(vary1)==0)         {
      vary1 <- n
      vary_name1 <- "n"
    } else {
      vary2 = n = t(n)
      vary_name2 = "n"
    }
  }

  if (is.numeric(MDRI)) {MDRI <- MDRI/365.25}
  if (is.numeric(BigT)) {BigT <- BigT/365.25}

  if (n=="out") {
    PrevR <- ((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)
    out2 <- PrevHR <- PrevH*PrevR
    out3 <- PrevHnR <-PrevH-PrevHR

#call DM_FirstOrderTerms instead of defined here.
    fot_PrevH <- ((((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)-FRR)/(((1-PrevH)^2)*(MDRI-FRR*BigT)))
    fot_PrevR <- (PrevH/((1-PrevH)*(MDRI-FRR*BigT)))
    fot_MDRI  <- ((FRR*PrevH-((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)*PrevH)/((1-PrevH)*((MDRI-FRR*BigT)^2)))
    fot_FRR   <- ((PrevH*(BigT*((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)-MDRI))/((1-PrevH)*((MDRI-FRR*BigT)^2)))

    out4 <- RSE_I_inf <- sqrt((fot_MDRI*RSE_MDRI*MDRI)^2+(fot_FRR*RSE_FRR*FRR)^2)/I

    out1 <- n <- ((fot_PrevH^2)*PrevH*(1-PrevH)*DE_H + (fot_PrevR^2)*(PrevR*(1-PrevR)*DE_R/(CR*PrevH)))/((RSE_I^2-RSE_I_inf^2)*I^2)
    #   alternative formula for n without requirement of previouscalculations (usable for uniroot)
    #   n <- (I^2*DE_H/(PrevH*(1-PrevH)) +
    #         DE_R*PrevH/((1-PrevH)^2*(MDRI-FRR*BigT)^2*CR)*
    #        (FRR+I*(1-PrevH)*(MDRI-FRR*BigT)/PrevH)*(1-(FRR+I*(1-PrevH)*(MDRI-FRR*BigT)/PrevH)))/
    #        ((RSE_I^2*I^2)-
    #        (((I*(1-PrevH)*(MDRI-FRR*BigT)*RSE_MDRI*MDRI)^2+
    #        ((PrevH*BigT*FRR+BigT*I*(1-PrevH)*(MDRI-FRR*BigT)-PrevH*MDRI)*FRR*RSE_FRR)^2)/
    #        (((1-PrevH)*(MDRI-FRR*BigT)^2)^2)))
    out5 <- RSE_PrevH <- sqrt(((PrevH*(1-PrevH))/n)*DE_H)/PrevH
    out6 <- RSE_PrevR <- sqrt(((PrevR*(1-PrevR))/n*CR*PrevH)*DE_R)/PrevR
    out_names <- c("sample.size","Prev.HIV&recent","Prev.HIV&nonrecent","RSE.I.inf.sample","RSE.PrevH", "RSE.PrevR")}


if(RSE_I=="out") {
    PrevR <- ((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)
    out2 <- PrevHR <- PrevH*PrevR
    out3 <- PrevHnR <-PrevH-PrevHR
###omit the function here defined below and input the ouptut of DM_FirstOrderTerms
    fot_PrevH <- ((((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)-FRR)/(((1-PrevH)^2)*(MDRI-FRR*BigT)))
    fot_PrevR <- (PrevH/((1-PrevH)*(MDRI-FRR*BigT)))
    fot_MDRI  <- ((FRR*PrevH-((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)*PrevH)/((1-PrevH)*((MDRI-FRR*BigT)^2)))
    fot_FRR   <- ((PrevH*(BigT*((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)-MDRI))/((1-PrevH)*((MDRI-FRR*BigT)^2)))

    out4 <- RSE_I_inf <- sqrt((fot_MDRI*RSE_MDRI*MDRI)^2+(fot_FRR*RSE_FRR*FRR)^2)/I

    out5 <- RSE_PrevH <- sqrt(((PrevH*(1-PrevH))/n)*DE_H)/PrevH
    out6 <- RSE_PrevR <- sqrt(((PrevR*(1-PrevR))/n*CR*PrevH)*DE_R)/PrevR

    out1 <- RSE_I <-sqrt(((fot_PrevH^2)*PrevH*(1-PrevH)*DE_H +
                            (fot_PrevR^2)*(PrevR*(1-PrevR)*DE_R/(CR*PrevH)))/(n*I^2)+RSE_I_inf^2)
    out_names <- c("RSE_I","Prev.HIV&recent","Prev.HIV&nonrecent","RSE.I.inf.sample","RSE.PrevH", "RSE.PrevR")}


  # if(PrevH=="out")    { }
  # if (CR=="out")      { }
  # if (MDRI=="out")    { }
  # if(RSE_MDRI=="out") { }
  # if (FRR=="out")     { }
  # if (RSE_FRR=="out") { }
  # if (BigT=="out")    { }
  # if (DE_H=="out")    { }
  # if(DE_R=="out")     { }

  if (length(test)==15) {
    VARY1<-vector(length=step)
    for (i in c(1:step)) {
      VARY1[i]<-paste(vary_name1,vary1[i,1])
    }
    VARY2<-vector(length=step)
    for (i in c(1:step)) {
      VARY2[i]<-paste(vary_name2,vary2[1,i])
    }
    VARY2 <- c("",VARY2)
    out1<-cbind(VARY1,out1)
    out1<-rbind(VARY2,out1)
    out2<-cbind(VARY1,out2)
    out2<-rbind(VARY2,out2)
    out3<-cbind(VARY1,out3)
    out3<-rbind(VARY2,out3)
    out4<-cbind(VARY1,out4)
    out4<-rbind(VARY2,out4)
    out5<-cbind(VARY1,out5)
    out5<-rbind(VARY2,out5)
    out6<-cbind(VARY1,out6)
    out6<-rbind(VARY2,out6)
  }

  if (length(test)==14) {
    VARY1<-vector(length=step)
    for (i in c(1:step)) {
      VARY1[i]<-paste(vary_name1,vary1[i,1])
    }
    out1<-cbind(VARY1,out1[,1])
    out2<-cbind(VARY1,out2[,1])
    out3<-cbind(VARY1,out3[,1])
    out4<-cbind(VARY1,out4[,1])
    out5<-cbind(VARY1,out5[,1])
    out6<-cbind(VARY1,out6[,1])
  }

  output <- list (out1, out2, out3, out4, out5,out6)
  names(output)<-out_names

  return(output)

}







#Examples of function use:
#####################################################################################################################
SSCprecision             ( I              =0.015,
                           RSE_I          =0.25,
                           PrevH          =0.2,
                           CR             =1,
                           MDRI           =200,
                           RSE_MDRI       =0.05,
                           FRR            =0.01,
                           RSE_FRR        =0.2,
                           BigT           = 730,
                           DE_H           = 1,
                           DE_R           = 1,
                           n              = "out",
                           step           = 5)
#####################################################################################################################
SSCprecision             ( I              =0.015,
                           RSE_I          =0.25,
                           PrevH          =c(0.15,0.25),
                           CR             =1,
                           MDRI           =200,
                           RSE_MDRI       =0.05,
                           FRR            =0.01,
                           RSE_FRR        =0.2,
                           BigT           = 710,
                           DE_H           = 1,
                           DE_R           = 1,
                           n              = "out",
                           step           = 5)


#####################################################################################################################
SSCprecision             ( I              =0.015,
                           RSE_I          ="out",
                           PrevH          =0.2,
                           CR             =0.7,
                           MDRI           =200,
                           RSE_MDRI       =0.05,
                           FRR            =0.01,
                           RSE_FRR        =0.2,
                           BigT           = 730,
                           DE_H           = 1,
                           DE_R           = 1,
                           n              = 3622,
                           step           = 5)



#####################################################################################################################
#there is some error here, not sure what...
SSCprecision             ( I              =0.015,
                           RSE_I          ="out",
                           PrevH          =0.2,
                           CR             =c(0.7,1),
                           MDRI           =200,
                           RSE_MDRI       =0.05,
                           FRR            =0.01,
                           RSE_FRR        =0.2,
                           BigT           = 730,
                           DE_H           = c(1,2),
                           DE_R           = 1,
                           n              = 3622,
                           step           = 5)




