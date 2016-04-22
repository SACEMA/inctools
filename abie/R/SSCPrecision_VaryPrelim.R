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
#
#
MDRI     <- 200
RSE_MDRI <- 0.05
FRR      <- 0.01
RSE_FRR  <- 0.2
BigT     <- 730
I        <- 0.015
RSE_I    <- "out"
PrevH    <- 0.1
CR       <- 1
DE_H     <- 1
DE_R     <- c(1,1.5)
n        <- c(5000,5500)
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


#Needs to be defined for function to work:
DM_FirstOrderTerms <- function (prevH, prevR, mdri, frr, bigt)   {
  fot_prevH <- (prevR-frr)/(((1-prevH)^2)*(mdri-frr*bigt)) #E.G. d(I)/d(P_H)
  fot_prevR <- prevH/((1-prevH)*(mdri-frr*bigt))
  fot_mdri  <- (frr*prevH-prevR*prevH)/((1-prevH)*((mdri-frr*bigt)^2))
  fot_frr   <- (prevH*(bigt*prevR-mdri))/((1-prevH)*((mdri-frr*bigt)^2))
  return (c(fot_prevH, fot_prevR, fot_mdri, fot_frr))
}



#' Sample size or precision calculation
#'
#' @param I Expected Incidence.
#' @param RSE_I Relative Standard Error of Incidence Estimate.
#' @param PrevH Prevalence of HIV.
#' @param CR Coverage rate: probability (0-1) of being tested for recency when positive for HIV.
#' @param MDRI mean duration of recent infection in days (vector/integer).
#' @param RSE_MDRI Relative standard error of MDRI (vector/integer).
#' @param FRR False recent rate (vector/integer).
#' @param RSE_FRR Relative standard error of FRR (vector/integer).
#' @param BigT post-infection time cut-off for true vs. false recency. Default is 730 days.
#' @param DE_H Design effect of HIV prevalence test (vector/integer).
#' @param DE_R Design effect of recency test (vector/integer).
#' @param n Sample Size: either a given hypothetical value, or to be determined by function, which is the default.
#' @param step number of steps between minimum I and maximum I in the calculation of a range of output.
#' @return Either sample size necessary for a given precision under a given set of testing characteristics and a hypothetical prevalence/incidence scenario, or precision under a particular sample size scenario, with a given hypothetical prevalence/incidence scenario.
#' @details
#'
#' Summarizes performance of a recent infection test (into a standard error of the incidence estimate), given estimated test properties (RSE of incidence) and the prevalence/incidence in a hypothetical context; or gives sample size necessary for a given level of estimator precision.
#' Returns: proportion of sample categorized as HIV positive and recently infected; proportion of sample categorized as HIV positive and non-recently infected; the relative standard error of the incidence estimator at infinite sample size, which is the component of variability explained solely by the assay characteristics; the relative standard error of the estimate of prevalence; the relative standard error of the estimate of proportion of HIV positive that are recent.
#'
#' @examples
#' SSCprecision(I = 0.015, RSE_I = 0.25, PrevH = 0.2, CR = 1,
#' MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2,
#' BigT = 730, DE_H = 1.1, DE_R = 1, n = "out")
#'
#' SSCprecision(I = c(0.015,0.02), RSE_I = 0.25, PrevH = c(0.10,0.20),
#' CR = 1, MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2,
#' BigT = 700, DE_H = 1, DE_R = 1, n = "out", step = 5)
#'
#' SSCprecision(I = 0.017, RSE_I = "out", PrevH = c(0.10,0.20),
#' CR = 1, MDRI = 211, RSE_MDRI = 0.05, FRR = 0.009, RSE_FRR = 0.2,
#' BigT = 720, n = 5000, step = 5)
#' @export
#'
#FOR THIS FUNCTION TO RUN, THE FUNCTION DM_FirstOrderTerms MUST BE INVOKED. IT EXISTS IN R SCRIPT recencyI_prev.R
SSCprecision <-     function ( I              ,
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
  var_list <- list(I=I, RSE_I=RSE_I, PrevH=PrevH, CR=CR, MDRI=MDRI, RSE_MDRI=RSE_MDRI, FRR=FRR, RSE_FRR=RSE_FRR, BigT=BigT, DE_H=DE_H, DE_R=DE_R, n=n, step=step)

  #CHECK TO MAKE SURE ONLY TWO VARIABLES ARE ALLOWED TO VARY
  max.list<-0
  for (i in 1:length(var_list)){
   if (length(var_list[[i]]) > 1) {max.list=max.list+1}
   if (max.list>2)  {stop("only a maximum of 2 variables are allowed to vary")}
  }

  if (sum(var_list=="out") > 1) {
    stop("only one of the variables RSE_I or n can be requested at a time")
  }

  if (length(var_list[1:8]) < 8) {
    stop("Not enough variables have been specified")
  }

for (i in 1:12) {
  if (length(var_list[[i]])>2 | length(var_list[[i]])<1)  {
    stop(paste("specifiy (only) min & max values for ",names(var_list)[i]),sep="")
    }
}

for (i in c(1,3,4,6:8)) {
  if (is.numeric(var_list[[1:length(var_list[i])]]) > 0 &
      (sum(var_list[[i]]<=1)!=length(var_list[[i]]) | sum(var_list[[i]]>=0)!=length(var_list[[i]]))  )
  {stop("Some input values are less than 0 or greater than 1")}
}

if(sum(RSE_I!="out") > 0) {
  if (is.numeric(var_list[[1:length(var_list[2])]]) > 0 &
      (sum(var_list[[2]]<=1)!=length(var_list[[2]]) | sum(var_list[[2]]>=0)!=length(var_list[[2]]))  )
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

  if (sum(BigT <= 182)){
    warning ("BigT is smaller than half a year")
      }
  if (sum(BigT < MDRI) > 0){
    stop ("MDRI cannot be greater than BigT")
  }
  if(sum(RSE_MDRI < 0.01) > 0){
    warning("RSE of estimated MDRI is less than 1%")
  }
  if(sum(FRR==0) > 0){
    warning("Zero estimated FRR")
  }
  if(sum(FRR > 0.10) > 0){
    warning("Estimated FRR is greater than 10%")
  }
  if(sum(RSE_FRR > 0.30) > 0){
    warning("RSE of estimated FRR is greater than 30%")
  }
  if(sum(RSE_FRR < 0.05) > 0){
    warning("RSE of estimated FRR is less than 5%")
  }
  if(sum(I > 0.20) > 0){
    warning(paste("Possible error in incidence input.", max(I)," seems exceptionally high",sep=""))
  }
  if(sum(n < 1000) > 0) {
    warning("Sample size is smaller than 1000")
  }


#Need error that reads: "ERROR: Test properties not consistent with test for recent infection"


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


  PrevR <- ((I*(1-PrevH)*(MDRI-FRR*BigT)) / PrevH + FRR)
  out2  <- PrevHR  <- round(PrevH*PrevR,  digits=5)     #Prev.HIV&recent
  out3  <- PrevHnR <- round(PrevH-PrevHR, digits=5)    #Prev.HIV&nonrecent

  fot<-DM_FirstOrderTerms(PrevH, PrevR, MDRI, FRR, BigT)
  #if the output of each term of DM_FirstOrderTerms is univariate, do one thing, otherwise, do another...
  if (sum(lengths(var_list) > 1)==0 | length(fot)==4) {
    fot_PrevH <- fot[1]
    fot_PrevR <- fot[2]
    fot_MDRI  <- fot[3]
    fot_FRR   <- fot[4]
  } else { #this needs to reflect what happens to the output fot matrix which depends on which variables are allowed to vary.
    if (length(I)>1 & sum(lengths(var_list))==14){
       fot_PrevH <- matrix(fot[1:(step*step)],nrow=step,ncol=step,byrow=F)
       fot_PrevR <- fot[(step*step+1)] #things change if and only if only Incidence is allowed to vary.
       fot_MDRI  <- matrix(fot[((step*step)+2):(((step*step)*2+1))],nrow=step,ncol=step)
       fot_FRR   <- matrix(fot[(((step*step)*2+2)):length(fot)],nrow=step,ncol=step)
    }else{ #here is the situation if only a variable besides incidence is allowed to vary, or any two parameters are allowed to vary.
    fot_PrevH <- matrix(fot[1:(step*step)],nrow=step,ncol=step,byrow=F)
    fot_PrevR <- matrix(fot[((step*step)*1+1):(((step*step)*2))],nrow=step,ncol=step)
    fot_MDRI  <- matrix(fot[((step*step)*2+1):(((step*step)*3))],nrow=step,ncol=step)
    fot_FRR   <- matrix(fot[((step*step)*3+1):(((step*step)*4))],nrow=step,ncol=step)
    }
  }


#IF SAMPLE SIZE n IS THE OUPUT VARIABLE (SO PRECISION/RSE_I IS FIXED)
  if (sum(n=="out")>0) {
    out4 <- RSE_I_inf_ss <- round(sqrt((fot_MDRI*RSE_MDRI*MDRI)^2+(fot_FRR*RSE_FRR*FRR)^2)/I, digits=5) #RSE.I.inf.sample
    out1 <- n <- ceiling(((fot_PrevH^2)*PrevH*(1-PrevH)*DE_H + (fot_PrevR^2)*(PrevR*(1-PrevR)*DE_R/(CR*PrevH)))/((RSE_I^2-RSE_I_inf_ss^2)*I^2))

    if(sum(((PrevH*(1-PrevH))/n)*DE_H<=0)>0){
      stop("no sample size will meet input constraints")
      }

    out5 <- RSE_PrevH <- round(sqrt(((PrevH*(1-PrevH))/n)*DE_H)/PrevH, digits=5)              #RSE.PrevH
    out6 <- RSE_PrevR <- round(sqrt(((PrevR*(1-PrevR))/n*CR*PrevH)*DE_R)/PrevR, digits=5)     #RSE.PrevR
    out_names <- c("sample.size","Prev.HIV.and.recent","Prev.HIV.and.nonrecent","RSE.I.inf.sample","RSE.PrevH", "RSE.PrevR")
    }


#IF precision IS THE OUPUT VARIABLE (SO sample size n IS FIXED)
if (sum(RSE_I=="out") > 0) {
    out4 <- RSE_I_inf_ss <- round(sqrt((fot_MDRI*RSE_MDRI*MDRI)^2+(fot_FRR*RSE_FRR*FRR)^2)/I, digits=5)

    out5 <- RSE_PrevH <- round(sqrt(((PrevH*(1-PrevH))/n)*DE_H)/PrevH, digits=5)
    out6 <- RSE_PrevR <- round(sqrt(((PrevR*(1-PrevR))/n*CR*PrevH)*DE_R)/PrevR, digits=5)

    out1 <- RSE_I <-round(sqrt(((fot_PrevH^2)*PrevH*(1-PrevH)*DE_H + (fot_PrevR^2)*(PrevR*(1-PrevR)*DE_R/(CR*PrevH)))/(n*I^2)+RSE_I_inf_ss^2), digits=5)
    out_names <- c("RSE_I","Prev.HIV.and.recent","Prev.HIV.and.nonrecent","RSE.I.inf.sample","RSE.PrevH", "RSE.PrevR")
    }
#
#   if (sum(RSE_I > 0.50) >0) {
#     warning("Implied RSE of incidence is greater than 50%")
#     }


  if (sum(lengths(var_list))==15) { #if two variables are to vary
    variable.1 <- vector(length=step)
    variable.2 <- vector(length=step)
    for (i in 1:step) {
      variable.1[i] <- paste(vary_name1, "=", vary1[i,1])
      variable.2[i] <- paste(vary_name2, "=", vary2[1,i])
    }

    if(length(out1)>1){
      out1 <- data.frame(out1)
      row.names(out1) <-variable.1
      colnames(out1)  <-variable.2
    }
    if(length(out2)>1){
      out2 <- data.frame(out2)
      row.names(out2) <-variable.1
      colnames(out2)  <-variable.2
    }
    if(length(out3)>1){
      out3 <- data.frame(out3)
      row.names(out3) <-variable.1
      colnames(out3)  <-variable.2
    }
    if(length(out4)>1){
      out4 <- data.frame(out4)
      row.names(out4) <-variable.1
      colnames(out4)  <-variable.2
    }
    if(length(out5)>1){
      out5 <- data.frame(out5)
      row.names(out5) <-variable.1
      colnames(out5)  <-variable.2
    }
    if(length(out6)>1){
      out6 <- data.frame(out6)
      row.names(out6) <-variable.1
      colnames(out6)  <-variable.2
    }


    # out1<-data.frame(out1)
    # out2<-data.frame(out2)
    # out3<-data.frame(out3)
    # out4<-data.frame(out4)
    # out5<-data.frame(out5)
    # out6<-data.frame(out6)
    #
    # row.names(out1) <-variable.2
    # colnames(out1)  <-variable.1
    # row.names(out2) <-variable.2
    # colnames(out2)  <-variable.1
    # row.names(out3) <-variable.2
    # colnames(out3)  <-variable.1
    # row.names(out4) <-variable.2
    # colnames(out4)  <-variable.1
    # row.names(out5) <-variable.2
    # colnames(out5)  <-variable.1
    # row.names(out6) <-variable.2
    # colnames(out6)  <-variable.1
  }


  if (sum(lengths(var_list))==14) { #if just one variable is allowed to vary
    variable.1 <- vector(length=step)
    for (i in c(1:step)) {
      variable.1[i] <- paste(vary_name1,"=", vary1[i,1])
    }
    if(length(out1)>1){
      out1 <- data.frame(variable.1,out1[,1])
      names(out1) <- c(vary_name1,out_names[1])
    }
    if(length(out2)>1){
      out2 <- data.frame(variable.1,out2[,1])
      names(out2) <- c(vary_name1,out_names[2])
    }
    if(length(out3)>1){
      out3 <- data.frame(variable.1,out3[,1])
      names(out3) <- c(vary_name1,out_names[3])
    }

    if(length(out4)>1){
      out4 <- data.frame(variable.1,out4[,1])
      names(out4) <- c(vary_name1,out_names[4])
    }
    if(length(out5)>1){
      out5 <- data.frame(variable.1,out5[,1])
      names(out5) <- c(vary_name1,out_names[5])
    }
    if(length(out6)>1){
      out6 <- data.frame(variable.1,out6[,1])
      names(out6) <- c(vary_name1,out_names[6])
    }
  }

  output <- list (out1, out2, out3, out4, out5, out6)
  names(output) <- out_names

  return(output)
  }







#Examples of function use:
#####################################################################################################################
#default of spreadsheet ABIE_v3_Sample_Size_Calculator.xlsx

SSCprecision             ( I              =0.015,
                           RSE_I          =0.25,
                           PrevH          =0.20,
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
                           PrevH          =0.20,
                           CR             =1,
                           MDRI           =200,
                           RSE_MDRI       =0.05,
                           FRR            =0.01,
                           RSE_FRR        =0.2,
                           BigT           = c(530,730),
                           DE_H           = 1,
                           DE_R           = 1,
                           n              = "out",
                           step           = 5)


SSCprecision             ( I              =c(0.015,.02),
                           RSE_I          =0.25,
                           PrevH          =0.20,
                           CR             =1,
                           MDRI           =200,
                           RSE_MDRI       =0.05,
                           FRR            =0.01,
                           RSE_FRR        =0.2,
                           BigT           = c(530,730),
                           DE_H           = 1,
                           DE_R           = 1,
                           n              = "out",
                           step           = 5)



#doesn't work when FRR goes above 3.7% for these values
# SSCprecision             ( I              =0.015,
#                            RSE_I          =0.25,
#                            PrevH          =0.20,
#                            CR             =1,
#                            MDRI           =200,
#                            RSE_MDRI       =0.05,
#                            FRR            =0.039,
#                            RSE_FRR        =0.2,
#                            BigT           = 730,
#                            DE_H           = 1,
#                            DE_R           = 1,
#                            n              = "out",
#                            step           = 5)

#####################################################################################################################
#default of spreadsheet ABIE_v3_Test_Performance_Calculator
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
                           n              = c(5000,5500),
                           step           = 5)
#this does not give same value of RSE_I as spreadhsheet for these values, close, but not totally
#this is not my (Avery's) error: but was present in sheet Petra gave me.

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
SSCprecision             ( I              =0.015,
                           RSE_I          ="out",
                           PrevH          =c(0.2,0.22),
                           CR             =1,
                           MDRI           =200,
                           RSE_MDRI       =0.05,
                           FRR            =0.01,
                           RSE_FRR        =0.2,
                           BigT           = 730,
                           DE_H           = c(1,1.1),
                           DE_R           = 1,
                           n              = 3622,
                           step           = 5)








#there is some error here, not sure what...
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




