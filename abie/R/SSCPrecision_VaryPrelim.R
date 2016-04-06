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

MDRI     <- 200
RSE_MDRI <- 0.05
FRR      <- 0.01
RSE_FRR  <- 0.2
BigT     <- 730
I        <- 0.015
RSE_I    <- 0.25
PrevH    <- 0.1
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
# and other variables of test, hypothetical data, and gives precision), and what was done in sheet
# ABIE_v3_Sample_Size_Calculator, which gives SS for a given precision level.




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
#' @param n Sample Size, either a given hypothetical value or to be determined by function (default)
#' @param step number of steps between minimum I and maximum I in the calculation of a range of output
#' @return Either sample size necessary for a given precision under a given set of testing characteristics and a hypothetical prevalence/incidence scenario, or precision under a particular sample size scenario, with a given hypothetical prevalence/incidence scenario.
#' @details
#'
#' Summarizes performance of a recent infection test (into a standard error of the incidence estimate), given estimated test properties and the prevalence/incidence in a hypothetical context, or gives sample size necessary for a given level of estimator precision.
#' Returns: proportion of sample categorized as HIV positive and recently infected; proportion of sample categorized as HIV positive and non-recently infected; the relative standard error of the incidence estimator at infinite sample size, which is the component of variability explained soley by the assay characteristics; the relative standard error of the estimate of prevalence; the relative standard error of the estimate of proportion of HIV positive that are recent.
#'
#' @examples
#' SSCprecision(I = 0.015, RSE_I = 0.25, PrevH = 0.2, CR = 1,
#' MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2,
#' BigT = 730, DE_H = 1, DE_R = 1, n = "out", step = 5)
#' @export
#'
#FOR THIS FUNCTION TO RUN, THE FUNCTION DM_FirstOrderTerms MUST BE INVOKED. IT EXISTS IN R SCRIPT recencyI_prev.R
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
#still need to make sure the function only handles a single "out" and not two, or none
#CHECK TO MAKE SURE ONLY TWO VARIABLES ARE ALLOWED TO VARY
  var_list <- list(I=I, RSE_I=RSE_I, PrevH=PrevH, CR=CR, MDRI=MDRI, RSE_MDRI=RSE_MDRI, FRR=FRR, RSE_FRR=RSE_FRR, BigT=BigT, DE_H=DE_H, DE_R=DE_R, n=n, step=step)

  max.list<-0
  for(i in 1:length(var_list)){
  if (length(var_list[[i]]) > 1) {max.list=max.list+1}
  if (max.list>2)  {stop("only a maximum of 2 variables are allowed to vary")}
  }

  if (sum(var_list=="out") > 1) {
    stop("only one of the variables RSE_I or n can be requested at a time")
  }
  if (length(var_list[1:8]) < 8) {
    stop("Not enough variables have been specified")
  }

for(i in 1:12){
  if (length(var_list[[i]])>2 | length(var_list[[i]])<1)  {stop(paste("specifiy (only) min & max values for ",names(var_list)[i]),sep="")}
}


for(i in c(1:8,10,11)){
  if (is.numeric(var_list[[1:length(var_list[i])]]) > 0 &
      (sum(var_list[[i]]<=1)!=length(var_list[[i]]) | sum(var_list[[i]]>=0)!=length(var_list[[i]]))  )
  {stop("Some input values are less than 0 or greater than 1")}
}
#above code does what below code does, only in 1 line. Which should we keep?
  # if (is.numeric(I)>0)        {stopifnot (I<=1        & I>=0)}
  # if (is.numeric(RSE_I)>0)    {stopifnot (RSE_I<=1    & RSE_I>=0)}
  # if (is.numeric(PrevH)>0)    {stopifnot (PrevH<=1    & PrevH>=0)}
  # if (is.numeric(CR)>0)       {stopifnot (CR<=1       & CR>=0)}
  # if (is.numeric(MDRI)>0)     {stopifnot (MDRI>=0)}
  # if (is.numeric(RSE_MDRI)>0) {stopifnot (RSE_MDRI<=1 & RSE_MDRI>=0)}
  # if (is.numeric(FRR)>0)      {stopifnot (FRR<=1      & FRR>=0)}
  # if (is.numeric(RSE_FRR)>0)  {stopifnot (RSE_FRR<=1  & RSE_FRR>=0)}
  # if (is.numeric(DE_H)>0)     {stopifnot (DE_H>=1)}
  # if (is.numeric(DE_R)>0)     {stopifnot (DE_R>=1)}
  # if (is.numeric(n)>0)        {stopifnot (n>100)}


#WARNING FOR SIZE OF ASSAY TIME CUTPOINT
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

######### THIS WHOLE SECTION HERE IS FOR IF ONE OR MORE OF THE VARIABLES IS ALLOWED TO VARY#######
#CREATES TWO NULL VALUES, I BELIVE FOR WHICH VARIABLES ARE TO VARY
  vary1 <- NULL
  vary2 <- NULL
#NOW FOR EACH VARIABLE IN LIST, IF IT'S ALLOWED TO VARY, THEN MAKE A MATRIX OF VALUES FOR THE STEP,
#FROM BEGINNING TO END, STEP NUMBER OF COLUMNS. IF THE VARIABLE IS THE FIRST ONE IN THIS LIST TO VARY
#THEN WE'RE CALLING THE VARIABLE IT'S SAME NAME, AS A MATRIX, AND IF IT'S THE SECOND, WE'RE CALLING
#THAT VARIABLE AS IT'S SAME NAME, BUT AS A TRANSPOSED MATRIX, SO NOW THE ROWS ARE THE STEPPING VALUES
#there's got to be a way to make this more efficient, like if two have been met, STOP...
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
    #IF THIS VARIABLE IS TO VARY, MAKE A MATRIX OF DIM (STEP*STEP) WHERE EACH COLUMN IS THE SEQUENCE OF VALUES
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
  #IF THIS VARIABLE IS TO VARY, MAKE A MATRIX OF DIM (STEP*STEP) WHERE EACH COLUMN IS THE SEQUENCE OF VALUES
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
######### END SECTION FOR IF ONE OR MORE OF THE VARIABLES IS ALLOWED TO VARY#######

#NOW MAKE VARIABLES IN DAYS TO BE IN UNITS YEARS
  if (is.numeric(MDRI)) {MDRI <- MDRI/365.25}
  if (is.numeric(BigT)) {BigT <- BigT/365.25}


  t


  #need to put ifelse() in here to deal with matrix fot output vs. scalar fot output
  fot<-DM_FirstOrderTerms(PrevH, PrevR, MDRI, FRR, BigT)
  #if the output of each term of DM_FirstOrderTerms is univariate, do one thing, otherwise, do another...
  if(length(var_list) == 13){
    fot_PrevH <- fot[1]
    fot_PrevR <- fot[2]
    fot_MDRI  <- fot[3]
    fot_FRR   <- fot[4]
  }else{
    fot_PrevH <- matrix(fot[1:(step*step)],nrow=step,ncol=step,byrow=F)
    fot_PrevR <- matrix(fot[((step*step)*1+1):(((step*step)*2))],nrow=step,ncol=step)
    fot_MDRI  <- matrix(fot[((step*step)*2+1):(((step*step)*3))],nrow=step,ncol=step)
    fot_FRR   <- matrix(fot[((step*step)*3+1):(((step*step)*4))],nrow=step,ncol=step)
  }

#IF SAMPLE SIZE n IS THE OUPUT VARIABLE (SO PRECISION/RSE_I IS FIXED)
  if (n=="out") {
    PrevR <- ((I*(1-PrevH)*(MDRI-FRR*BigT)) / PrevH + FRR)
    out2 <- PrevHR <- PrevH*PrevR     #Prev.HIV&recent
    out3 <- PrevHnR <-PrevH-PrevHR    #Prev.HIV&nonrecent
    out4 <- RSE_I_inf_ss <- sqrt((fot_MDRI*RSE_MDRI*MDRI)^2+(fot_FRR*RSE_FRR*FRR)^2)/I #RSE.I.inf.sample
    out1 <- n <- ((fot_PrevH^2)*PrevH*(1-PrevH)*DE_H + (fot_PrevR^2)*(PrevR*(1-PrevR)*DE_R/(CR*PrevH)))/((RSE_I^2-RSE_I_inf_ss^2)*I^2)
    out5 <- RSE_PrevH <- sqrt(((PrevH*(1-PrevH))/n)*DE_H)/PrevH              #RSE.PrevH
    out6 <- RSE_PrevR <- sqrt(((PrevR*(1-PrevR))/n*CR*PrevH)*DE_R)/PrevR     #RSE.PrevR
    out_names <- c("sample.size","Prev.HIV&recent","Prev.HIV&nonrecent","RSE.I.inf.sample","RSE.PrevH", "RSE.PrevR")
    }



#IF SAMPLE SIZE n IS THE OUPUT VARIABLE (SO PRECISION/RSE_I IS FIXED)
if(RSE_I=="out") {
    PrevR <- ((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)
    out2 <- PrevHR <- PrevH*PrevR
    out3 <- PrevHnR <-PrevH-PrevHR
###omit the function here defined below and input the ouptut of DM_FirstOrderTerms
    fot_PrevH <- ((((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)-FRR)/(((1-PrevH)^2)*(MDRI-FRR*BigT)))
    fot_PrevR <- (PrevH/((1-PrevH)*(MDRI-FRR*BigT)))
    fot_MDRI  <- ((FRR*PrevH-((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)*PrevH)/((1-PrevH)*((MDRI-FRR*BigT)^2)))
    fot_FRR   <- ((PrevH*(BigT*((I*(1-PrevH)*(MDRI-FRR*BigT))/PrevH + FRR)-MDRI))/((1-PrevH)*((MDRI-FRR*BigT)^2)))

    out4 <- RSE_I_inf_ss <- sqrt((fot_MDRI*RSE_MDRI*MDRI)^2+(fot_FRR*RSE_FRR*FRR)^2)/I

    out5 <- RSE_PrevH <- sqrt(((PrevH*(1-PrevH))/n)*DE_H)/PrevH
    out6 <- RSE_PrevR <- sqrt(((PrevR*(1-PrevR))/n*CR*PrevH)*DE_R)/PrevR

    out1 <- RSE_I <-sqrt(((fot_PrevH^2)*PrevH*(1-PrevH)*DE_H +
                            (fot_PrevR^2)*(PrevR*(1-PrevR)*DE_R/(CR*PrevH)))/(n*I^2)+RSE_I_inf_ss^2)
    out_names <- c("RSE_I","Prev.HIV&recent","Prev.HIV&nonrecent","RSE.I.inf.sample","RSE.PrevH", "RSE.PrevR")
    }



  if (length(var_list)==15) {
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

  if (length(var_list)==14) {
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
# SSCprecision             ( I              =0.015,
#                            RSE_I          ="out",
#                            PrevH          =0.2,
#                            CR             =c(0.7,1),
#                            MDRI           =200,
#                            RSE_MDRI       =0.05,
#                            FRR            =0.01,
#                            RSE_FRR        =0.2,
#                            BigT           = 730,
#                            DE_H           = c(1,2),
#                            DE_R           = 1,
#                            n              = 3622,
#                            step           = 5)




