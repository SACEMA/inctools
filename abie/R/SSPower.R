#This script is a function that takes either a sample size and gives power for a test of two incidences,
#or takes power and gives the required sample size necessary

#Needs to be defined for function to work:
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


#Test values:
#I<-c(0.05,0.07)
I<-NA
Power="out"
MDRI     <- 200
RSE_MDRI <- 0.05
FRR      <- 0.01
RSE_FRR  <- 0.2
BigT     <- 730

PrevH    <- c(0.2,0.15)
PrevR    <- c(0.18,0.13) #could be vector
CR       <- 1
DE_H     <- 1
DE_R     <- 1
N        <- c(5000,5000)
step <- 5
BMest="same.test"

#add N_testR to calculation...

SSPower <- function (I=NA, PrevH, PrevR,
                      I1, I2, N, alpha=0.05, Power="out", DE_H=1, DE_R=1,
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

  # if(BS_Count<=0){
  #   stop("Bootstrap samples count must be positive integer")
  # }

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

  #these next two lines were just added...
  if (length(PrevH)==1) {PrevH <- rep(PrevH, times=no_s)} else {PrevH=PrevH}
  if (length(PrevR)==1) {PrevR <- rep(PrevR, times=no_s)} else {PrevR=PrevR}
  if (is.na(I[1])) {I <- rep(I, times=no_s)}

  stopifnot(no_s==length(PrevR)    & no_s==length(PrevR) &
              no_s==length(MDRI)    & no_s==length(RSE_MDRI)  & no_s==length(FRR)   &
              no_s==length(RSE_FRR)  & length(BigT)==1)


  MDRI<-MDRI/365.25
  BigT<-BigT/365.25




  if(is.na(I[1])){
  I_Est <- I_EST(prevH=PrevH, prevR=PrevR, mdri=MDRI, frr=FRR, bigt=BigT)
  } else {I_Est <- I}

  # deltaI_Est_Mat<-matrix(ncol=no_s, nrow=no_s) #empty matrix of dim (# surveys*# surveys)
  # for (i in c(1:no_s)) {
  #   for (j in c(1:no_s)) {
  #     deltaI_Est_Mat[j,i]<-I_Est[i]-I_Est[j]
  #   }
  # }
  # #now each element in deltaI_Est_Mat is the difference in I for each survey
  # #the ij element is I[i]-I[j].
  #
  #
  # deltaI_Est_Vec <- as.vector(deltaI_Est_Mat)
  # #makes above matrix into a vector

deltaI_Est <- I_Est[1]-I_Est[2]


  # #this section defines the empirical variance of each object to be bootsrapped for the mvtrnorm function
  # if (Boot==TRUE) {
  #   BS_Var_PrevH <- (RSE_PrevH*PrevH)^2
  #   BS_Var_PrevR <- (RSE_PrevR*PrevR)^2
  #   BS_Var_MDRI  <- (MDRI*RSE_MDRI)^2
  #   BS_Var_FRR   <- (FRR*RSE_FRR)^2
  #
  #   DM_Var_PrevH <- rep(0, times=no_s)
  #   DM_Var_PrevR <- rep(0, times=no_s)
  #   DM_Var_MDRI  <- rep(0, times=no_s)
  #   DM_Var_FRR   <- rep(0, times=no_s)
  # } else {
    # BS_Var_PrevH <- rep(0, times=no_s)
    # BS_Var_PrevR <- rep(0, times=no_s)
    # BS_Var_MDRI  <- rep(0, times=no_s)
    # BS_Var_FRR   <- rep(0, times=no_s)
    #
    # DM_Var_PrevH <- (RSE_PrevH*PrevH)^2
    # DM_Var_PrevR <- (RSE_PrevR*PrevR)^2

DM_Var_PrevH<-PrevH*(1-PrevH)*DE_H/N #need to rename these two...
DM_Var_PrevR<-PrevR*(1-PrevR)*DE_R/N

    DM_Var_MDRI  <- (MDRI*RSE_MDRI)^2
    DM_Var_FRR   <- (FRR*RSE_FRR)^2
#  }


  # if (Boot==TRUE) {
  #
  #   DM_Var_I <- rep(0, times=no_s)
  #   DM_Var_deltaI <- rep(0, times=no_s^2)
  #
  #   I_BSMat <- matrix(nrow=BS_Count, ncol=no_s) #creates empty matrix of dim (#BS samples*#surveys)
  #   BS_RootEstMat <- matrix(nrow=BS_Count, ncol=no_s*4) #creates empty matrix of dim (#BS samples-by-#surveys*4)
  #   BS_Var_I <- vector(length=no_s) #vector of length # of surveys
  #
  #   #loop that, for each survey
  #   for (i in 1:no_s) {
  #     BS_RootEstMat [,(i*4-3):(i*4)] <- BS_SURVEY_ESTS  (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i], frr=FRR[i], bs_count=BS_Count,
  #                                                        bs_var_prevH=BS_Var_PrevH[i], bs_var_prevR=BS_Var_PrevR[i],
  #                                                        bs_var_mdri=BS_Var_MDRI[i], bs_var_frr=BS_Var_FRR[i],
  #                                                        covar_HR=Covar_HR[i])
  #   }
  #   #returns bootstraps of prevH, prevR, mdri, frr as columns, one for each survey, so if survey#=3, then
  #   #first four columns are for prevH, prevR, mdri, frr as columns, then next four are for those varaibles from second survey, and so on...
  #
  #   #I THINK THIS BELOW SECTION SHOULD HAVE AN INDEX LOOP RIGHT?? AS IT STANDS THERE'S NOTHING...
  #   for(i in 1:no_s){
  #     if ((BMest=="same.test"| BMest=="FRR.indep")  ) {
  #       BS_RootEstMat[,(i*4-1)] <- BS_RootEstMat[,3]
  #     }#above is saying, if the test is the same, or FRR.indep (meaning mdri still same) then each column in BS matrix
  #     #corresponding to mdri will equal the first BS sample of that variable
  #     if (BMest=="same.test" ) {
  #       BS_RootEstMat[,(i*4)] <- BS_RootEstMat[,4]
  #     } #similarly if it's same test then FRR will be recorded as same in there.
  #   }
  #
  #
  #   for (i in 1:no_s) {
  #     I_BSVec <- I_EST(prevH=BS_RootEstMat[,(i*4-3)], prevR=BS_RootEstMat[,(i*4-2)],
  #                      mdri=BS_RootEstMat[,(i*4-1)],  frr=BS_RootEstMat[,(i*4)], bigt=BigT)
  #     I_BSMat[,i] <-I_BSVec
  #     BS_Var_I[i] <- var(I_BSMat[,i])
  #   }
  #   #now compute BS I and variance of I from BS samples
  #
  #
  #   #matrix of bootstrapped differences between surveys
  #   deltaI_BSMat <- matrix(nrow=BS_Count, ncol=(no_s^2))#empty matrix of dim (BS samples by number surveys squared)
  #   for (i in c(1:no_s)){
  #     for (j in c(1:no_s)) {
  #       deltaI_BSMat[,(i*no_s-(no_s-j))] <- I_BSMat[,i]-I_BSMat[,j]
  #       I_BSMat[,i] <- sort(I_BSMat[,i], decreasing=FALSE)
  #     }
  #   }
  #
  #   #makes vector of variances of difference in BS estimates of I from different surveys
  #   BS_Var_deltaI <- vector(length=no_s^2)
  #   for (i in c(1:(no_s^2))) {
  #     deltaI_BSMat[,i] <- sort(deltaI_BSMat[,i], decreasing=FALSE)
  #     BS_Var_deltaI[i] <- var(deltaI_BSMat[,i])
  #   }
  #
  #   #tracks back to line: 'if (Boot==TRUE) {'
  #
  # }
  # else { # if Boot is FALSE
    # BS_Var_I=rep(0, times=no_s) #so if no BS-ing, make the BS variables null
    # BS_Var_deltaI=rep(0, times=no_s^2)


    #next few lines make delta method matrix
    fot_Mat  <- matrix (nrow=no_s, ncol=4)
    DM_Var_I <- vector(length=no_s)
    for (i in c(1:no_s)) {
      fot_Mat[i,] <- DM_FirstOrderTerms (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i], frr=FRR[i], bigt=BigT)
      #DM_FirstOrderTerms gives fot_prevH, fot_prevR, fot_mdri, fot_frr, so FOT for each prev. survey input
      #and each row of fot_Mat has the FOT for prev. survey i.
      DM_Var_I[i] <- (fot_Mat[i,1]^2)*DM_Var_PrevH[i] + (fot_Mat[i,2]^2)*DM_Var_PrevR[i] +
        (fot_Mat[i,3]^2)*DM_Var_MDRI[i]  + (fot_Mat[i,4]^2)*DM_Var_FRR[i]
    }

    # DM_Var_deltaI <- vector(length=no_s^2)
    # for (i in c(1:no_s)) {
    #   for (j in c(1:no_s)) {
        # DM_Var_deltaI[i*no_s-(no_s-j)] <- DM_VAR_deltaI (BMest=BMest, fot_prevH1=fot_Mat[i,1], fot_prevH2=fot_Mat[j,1],
        #                                                  fot_prevR1=fot_Mat[i,2],   fot_prevR2=fot_Mat[j,2],
        #                                                  fot_mdri1=fot_Mat[i,3],    fot_mdri2=fot_Mat[j,3],
        #                                                  fot_frr1=fot_Mat[i,4],     fot_frr2=fot_Mat[j,4],
        #                                                  dm_var_prevH1=DM_Var_PrevH[i], dm_var_prevH2=DM_Var_PrevH[j],
        #                                                  dm_var_prevR1=DM_Var_PrevR[i], dm_var_prevR2=DM_Var_PrevR[j],
        #                                                  dm_var_mdri1=DM_Var_MDRI[i],   dm_var_mdri2=DM_Var_MDRI[j],
        #                                                  dm_var_frr1=DM_Var_FRR[i],     dm_var_frr2=DM_Var_FRR[j])
    #   }
    # }

    DM_Var_deltaI <- DM_VAR_deltaI (BMest=BMest, fot_prevH1=fot_Mat[1,1], fot_prevH2=fot_Mat[2,1],
                                                     fot_prevR1=fot_Mat[1,2],   fot_prevR2=fot_Mat[2,2],
                                                     fot_mdri1=fot_Mat[1,3],    fot_mdri2=fot_Mat[2,3],
                                                     fot_frr1=fot_Mat[1,4],     fot_frr2=fot_Mat[2,4],
                                                     dm_var_prevH1=DM_Var_PrevH[1], dm_var_prevH2=DM_Var_PrevH[2],
                                                     dm_var_prevR1=DM_Var_PrevR[1], dm_var_prevR2=DM_Var_PrevR[2],
                                                     dm_var_mdri1=DM_Var_MDRI[1],   dm_var_mdri2=DM_Var_MDRI[2],
                                                     dm_var_frr1=DM_Var_FRR[1],     dm_var_frr2=DM_Var_FRR[2])


    DM_Var_deltaI.infSS <- DM_VAR_deltaI.infSS (BMest=BMest,
                                    fot_mdri1=fot_Mat[1,3],    fot_mdri2=fot_Mat[2,3],
                                    fot_frr1=fot_Mat[1,4],     fot_frr2=fot_Mat[2,4],
                                    dm_var_mdri1=DM_Var_MDRI[1],   dm_var_mdri2=DM_Var_MDRI[2],
                                    dm_var_frr1=DM_Var_FRR[1],     dm_var_frr2=DM_Var_FRR[2])


  DM_SD_I      <- sqrt(DM_Var_I)
  DM_SD_deltaI <- sqrt(DM_Var_deltaI)
  DM_SD_deltaI.infSS <- sqrt(DM_Var_deltaI.infSS)
  Var_I        <-  DM_Var_I
  RSE_I        <- sqrt(Var_I)/I_Est
  SD_I         <- sqrt(Var_I)
  Var_deltaI   <- DM_Var_deltaI #+  BS_Var_deltaI
  RSE_deltaI   <- sqrt(Var_deltaI)/abs(deltaI_Est) #sqrt(Var_deltaI)/abs(deltaI_Est_Vec)
  RSE_deltaI.infSS   <- sqrt(DM_Var_deltaI.infSS)/abs(deltaI_Est)
  SD_deltaI    <- sqrt(Var_deltaI)

  #this function takes bootstrap matrix, SD via delta method, and I estimates, and returns the spread version of those
  #IM REWRITING THIS SO THAT IT ONLY TAKES FULL BS OR DELTA-METHOD
  CI_BSandDM <- function (BSMat, DM_SD, Est) {
    if (sum(DM_Var_I)>0)  {
      for (i in c(1:length(Est))) {
        CI_Mat[i,1] <- qnorm(alpha/2, mean=Est[i], sd=DM_SD[i])
        CI_Mat[i,2] <- qnorm(1-alpha/2, mean=Est[i], sd=DM_SD[i])
      }
    }
    # else {
    #   for (i in c(1:ncol(BSMat))) {
    #     CI_Mat[i,1] <- quantile(BSMat[,i], alpha/2)
    #     CI_Mat[i,2] <- quantile(BSMat[,i], 1-alpha/2)
    #   }
    # }
    return (CI_Mat)
  }

  CI_Mat <- matrix(nrow=no_s, ncol=2)
  CI_I_Mat <- CI_BSandDM (BSMat=I_BSMat, DM_SD=DM_SD_I, Est=I_Est)

  deltaI_CI<-NULL
  deltaI_CI[1] <- qnorm(alpha/2, mean=deltaI_Est, sd=DM_SD_deltaI)
  deltaI_CI[2] <- qnorm(1-alpha/2, mean=deltaI_Est, sd=DM_SD_deltaI)



  #p_value <- pnorm((-abs(deltaI_Est_Vec)/SD_deltaI), mean=0, sd=1)*2
  Power <- 1-pnorm(q=qnorm(1-alpha/2), mean=deltaI_Est, sd=1/RSE_deltaI)
  Power.infSS <-1-pnorm(q=qnorm(1-alpha/2), mean=deltaI_Est, sd=1/RSE_deltaI.infSS)
  # x <- seq(-4, 4, length=100)
  # hx <- dnorm(x)
  # plot(x,hx,type="l", lty=2)
  # curve(dnorm(x,mean=abs(deltaI_Est), sd=DM_SD_deltaI),col="red",add=T)

  # survey_no <- vector(length=no_s^2)
  # out_I_Est <- vector(length=no_s^2)
  # out_RSE_I <- vector(length=no_s^2)
  # out_CI_I_lo  <- vector(length=no_s^2)
  # out_CI_I_up <- vector(length=no_s^2)
  # delta_code <- vector(length=no_s^2)
  # for (i in c(1:no_s)) {
  #   survey_no  [(i*no_s-(no_s-1)):(i*no_s)] <- c(i, rep("", times=(no_s-1)))
  #   out_I_Est  [(i*no_s-(no_s-1)):(i*no_s)] <- c(round(I_Est[i], digit=5), rep("", times=(no_s-1)))
  #   out_RSE_I  [(i*no_s-(no_s-1)):(i*no_s)] <- c(round(RSE_I[i], digit=5), rep("", times=(no_s-1)))
  #   out_CI_I_lo[(i*no_s-(no_s-1)):(i*no_s)] <- c(round(CI_I_Mat[i,1], digit=5), rep("", times=(no_s-1)))
  #   out_CI_I_up[(i*no_s-(no_s-1)):(i*no_s)] <- c(round(CI_I_Mat[i,2], digit=5), rep("", times=(no_s-1)))
  #   for (j in c(1:no_s)) {
  #     delta_code [(i*no_s-(no_s-j))] <- paste(i, j, sep=" vs ")
  #   }
  # }

  # if(sum(out_RSE_I>0.25)){warning("RSE of incidence estimator greater than 25%")}
  #
  #
  # for (i in c(1:no_s)) {
  #   deltaI_Est_Vec[(i*no_s-(no_s-i))]<-NA
  #   RSE_deltaI[(i*no_s-(no_s-i))]<-NA
  #   p_value[(i*no_s-(no_s-i))]<-NA
  #   for (j in c(1:2)) {
  #     CI_deltaI_Mat[(i*no_s-(no_s-i)),j] <- NA
  #   }
  # }
  #
  #
  #
  # for(i in 1:length(out_CI_I_lo)){
  #   if((out_CI_I_lo[i]!="" & out_CI_I_lo[i]<0) | ( (out_CI_I_lo[i]!="" & out_CI_I_lo[i]>1)))
  #   {warning("CI out of [0,1] bounds"); break}
  # }
  # for(i in 1:length(out_CI_I_up)){
  #   if((out_CI_I_up[i]!="" & out_CI_I_up[i]<0) | ( (out_CI_I_up[i]!="" & out_CI_I_up[i]>1)))
  #   {warning("CI out of [0,1] bounds"); break}
  # }
  #
  #
  # out_deltaI_Est <- round(deltaI_Est, digit=5)#round(deltaI_Est_Vec, digit=5)
  # out_RSE_deltaI <- round(RSE_deltaI, digit=5)
  # out_p_value <- round(p_value,5)
  # out_CI_deltaI_Mat <- round(CI_deltaI_Mat, digit=5)
  #
  # out_p_value<-ifelse(out_p_value<0.001,"<0.0001",out_p_value)
  #
  # if (length(I_Est)==1) {
  #   output <- data.frame ("Incidence"=out_I_Est,
  #                         "CI lo"=out_CI_I_lo,
  #                         "CI up"=out_CI_I_up,
  #                         "RSE"=out_RSE_I) } else {
  #                           output <- data.frame ("survey"=survey_no,
  #                                                 "Incidence"=out_I_Est,
  #                                                 "CI lo"=out_CI_I_lo,
  #                                                 "CI up"=out_CI_I_up,
  #                                                 "RSE"=out_RSE_I,
  #                                                 "compare"=delta_code,
  #                                                 "Diff"=out_deltaI_Est,
  #                                                 "CI Diff lo"=out_CI_deltaI_Mat[,1],
  #                                                 "CI Diff up"=out_CI_deltaI_Mat[,2],
  #                                                 "RSE Diff"=out_RSE_deltaI,
  #                                                 "p-value"=out_p_value) }

if(Power="out"){
  output<- list(data.frame(Survey=c(1,2), I_EST=round(I_Est,3), RSE_I=round(RSE_I,3), CI.low=round(CI_I_Mat[,1],3), CI.up=round(CI_I_Mat[,2],3)),
              data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=round(RSE_deltaI.infSS,3), Power=round(Power,3), Power.infSS=round(Power.infSS,3)))
}

if(n="out"){
#fill in later
}

  return (output)
}








