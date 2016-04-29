## ---- echo=FALSE---------------------------------------------------------
SSPower <- function (I1, I2, PrevH1, PrevH2, n1 = "both", n2 = "both", alpha = 0.05, Power = 0.80, SS = "out", CR = 1, DE_H = 1, DE_R = 1, BMest = "same.test", MDRI, RSE_MDRI, FRR, RSE_FRR, BigT = 730){

############ Begin warning messages ################
   stopifnot (PrevH1<=1     & PrevH1>=0)
   stopifnot (PrevH2<=1     & PrevH2>=0)
   stopifnot (MDRI>=0)
   stopifnot (RSE_MDRI<=1  & RSE_MDRI>=0)
   stopifnot (FRR<=1       & FRR>=0)
   stopifnot (RSE_FRR<=1   & RSE_FRR>=0)

  if(sum(BMest==c("same.test", "FRR.indep", "MDRI.FRR.indep"))==0){
    stop("BMest option must be same.test, FRR.indep, or MDRI.FRR.indep")
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

  if (sum(BigT<=182)>0) {warning ("BigT is smaller than half a year")}


  no_s <- 2 #dimension of inputs (number of surveys)
  if (length(MDRI)==1)     {MDRI <- rep(MDRI, times=no_s)}         else {MDRI=MDRI}
  if (length(FRR)==1)      {FRR  <- rep(FRR, times=no_s)}          else {FRR=FRR}
  if (length(RSE_MDRI)==1) {RSE_MDRI <- rep(RSE_MDRI, times=no_s)} else {RSE_MDRI=RSE_MDRI}
  if (length(RSE_FRR)==1)  {RSE_FRR  <- rep(RSE_FRR, times=no_s)}  else {RSE_FRR=RSE_FRR}
  if (length(DE_H)==1)  {DE_H  <- rep(DE_H, times=no_s)}  else {DE_H=DE_H}
  if (length(DE_R)==1)  {DE_R  <- rep(DE_R, times=no_s)}  else {DE_R=DE_R}
  if (length(CR)==1)  {CR  <- rep(CR, times=no_s)}  else {CR=CR}
  if(BMest=="MDRI.FRR.indep"){if (length(BigT)==1)  {BigT  <- rep(BigT, times=no_s)}  else {BigT=BigT} }

  ############ End warning messages ################




  MDRI<-MDRI/365.25
  BigT<-BigT/365.25
  if(Power=="out") N <- c(n1,n2)
  PrevH <- c(PrevH1,PrevH2)
  I <- c(I1,I2)
  deltaI_Est <- I[1]-I[2]

  HIV.neg<- 1-PrevH
  PrevR <- I*(1-PrevH)*(MDRI-FRR*BigT)+(FRR*PrevH)

  DM_Var_MDRI  <- (MDRI*RSE_MDRI)^2
  DM_Var_FRR   <- (FRR*RSE_FRR)^2

  MDRI.CI <- 365.25*data.frame(CI.low=qnorm(alpha/2, mean=MDRI, sd=sqrt(DM_Var_MDRI)), CI.up=qnorm(1-alpha/2, mean=MDRI, sd=sqrt(DM_Var_MDRI)))
  FRR.CI <- data.frame(CI.low=qnorm(alpha/2, mean=FRR, sd=sqrt(DM_Var_FRR)),CI.up=qnorm(1-alpha/2, mean=FRR, sd=sqrt(DM_Var_FRR)))
  
if(Power=="out"){
  if(n1<1 | n2<1) stop("Sample size input must be a positive integer")

  Var_I <- I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))
                + (RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2
                +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2
                )
  RSE_I.infss <- sqrt(I^2*((RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2
                  +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2))/I
  RSE_I <- sqrt(Var_I)/I
  Var_delta_I <-
    if(BMest=="MDRI.FRR.indep"){Var_I[1]+Var_I[2]
    } else if(BMest=="FRR.indep") { sum(I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))))  +   ((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2  +    sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)
      } else if(BMest=="same.test"){ sum(I^2*((1/N)*(1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))))  +   ((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)   +   ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2) }

  RSE_deltaI <- sqrt(Var_delta_I)/abs(deltaI_Est)
  RSE_deltaI.infSS <-  if(BMest=="MDRI.FRR.indep"){sqrt(sum(((RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2 +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)*I^2))/deltaI_Est
  } else if(BMest=="FRR.indep") {   sqrt(((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2  + sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2))/deltaI_Est
  } else if(BMest=="same.test"){    sqrt(((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)   +   ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2))/deltaI_Est }



  CI.low <- qnorm(alpha/2, mean=I, sd=sqrt(Var_I))
  CI.up <- qnorm(1-alpha/2, mean=I, sd=sqrt(Var_I))

  deltaI_CI<-NULL
  deltaI_CI[1] <- qnorm(alpha/2, mean=deltaI_Est, sd=sqrt(Var_delta_I))
  deltaI_CI[2] <- qnorm(1-alpha/2, mean=deltaI_Est, sd=sqrt(Var_delta_I))

  ss.power <- 1-pnorm(q=qnorm(1-alpha/2), mean=1/RSE_deltaI, sd=1)
  Power.infSS <-1-pnorm(q=qnorm(1-alpha/2), mean=1/RSE_deltaI.infSS, sd=1)

  #if(ss.power<0.7){warning("Probability of correct inference less than 70%")}
  #this warning is in spreadsheets, but I don't like it. User expected to have an idea about power.


  if(BMest=="FRR.indep"){
  output <- list(Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)), CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                 Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                 Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), CI.low=round(MDRI.CI[1,1],3), CI.up=round(MDRI.CI[1,2],3)),
                 Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR,3), CI.low=round(FRR.CI[,1],3), CI.up=round(FRR.CI[,2],3)),
                 Implied.Subject.Counts=data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) )
} else
  if(BMest=="same.test") {output <- list(Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)), CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                      Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                      Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), CI.low=round(MDRI.CI[1,1],3), CI.up=round(MDRI.CI[1,2],3)),
                      Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR[1],3), CI.low=round(FRR.CI[1,1],3), CI.up=round(FRR.CI[1,2],3)),
                      Implied.Subject.Counts=data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) )
} else
  if(BMest=="MDRI.FRR.indep"){ output <- list(Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)), CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                                                  Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                                                  Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI*365.25,3), CI.low=round(MDRI.CI[,1],3), CI.up=round(MDRI.CI[,2],3)),
                                                  Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR,3), CI.low=round(FRR.CI[,1],3), CI.up=round(FRR.CI[,2],3)),
                                                  Implied.Subject.Counts=data.frame(Survey.1=c(HIV.negative=round(N[1]*HIV.neg[1]),HIV.positive=round(N[1]-N[1]*HIV.neg[1]),HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=round(N[2]*HIV.neg[2]),HIV.positive=round(N[2]-N[2]*HIV.neg[2]),HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) )
}
}

  else
    if(SS=="out" & n1=="out"){
      if(n2=="out"){stop("both n1 and n2 cannot be designated 'out'")}
    if (BMest=="same.test")
    {SS <- ceiling(I[1]^2*((1/PrevH[1])*(DE_H[1]/(1-PrevH[1])+(DE_R[1]/CR[1])*(PrevR[1]/PrevH[1])*(1-PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1]-FRR[1])^2)))/
                     ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)) )^2 -
                        ((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)  -
                        ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2) -
                        1/n2*I[2]^2*((1/PrevH[2])*(DE_H[2]/(1-PrevH[2])+(DE_R[2]/CR[2])*(PrevR[2]/PrevH[2])*(1-PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2]-FRR[2])^2)))
                     )
    )
    } else if(BMest=="FRR.indep")
    { SS <- ceiling(I[1]^2*((1/PrevH[1])*(DE_H[1]/(1-PrevH[1])+(DE_R[1]/CR[1])*(PrevR[1]/PrevH[1])*(1-PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1]-FRR[1])^2))) /
                      ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)))^2  -  ((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2 -
                         sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2) -
                         1/n2*I[2]^2*((1/PrevH[2])*(DE_H[2]/(1-PrevH[2])+(DE_R[2]/CR[2])*(PrevR[2]/PrevH[2])*(1-PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2]-FRR[2])^2)))
                      ))
    } else if(BMest=="MDRI.FRR.indep")
    { SS <- ceiling(I[1]^2*((1/PrevH[1])*(DE_H[1]/(1-PrevH[1])+(DE_R[1]/CR[1])*(PrevR[1]/PrevH[1])*(1-PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1]-FRR[1])^2))) /
                      ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)))^2  -  sum(I^2*(RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2) -
                         sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2 ) -
                         1/n2*I[2]^2*((1/PrevH[2])*(DE_H[2]/(1-PrevH[2])+(DE_R[2]/CR[2])*(PrevR[2]/PrevH[2])*(1-PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2]-FRR[2])^2)))
                      ))
    }

    if(SS<0){stop("No sample size will meet the given constraints")}


    #Now based on derived common SS, output implied summary statistics
    N<-c(SS,n2) #make derived necessary common sample size a vector
    if(sum(PrevR/PrevH*(CR*(N-N*HIV.neg))<10)>0 ) #formula is count of expected recent infections
      warning("Expected count of 'recent' infections is less than 10 for at least one survey")

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
    RSE_deltaI.infSS <-  if(BMest=="MDRI.FRR.indep"){ sqrt(sum(((RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2 +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)*I^2))/deltaI_Est
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
      output <- list(Minimum.SS=SS,
                     Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)), CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                     Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                     Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), CI.low=round(MDRI.CI[1,1],3), CI.up=round(MDRI.CI[1,2],3)),
                     Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR,3), CI.low=round(FRR.CI[,1],4), CI.up=round(FRR.CI[,2],4)),
                     Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
    } else
      if(BMest=="same.test") {output <- list(Minimum.SS=SS,
                                             Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)),CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                                             Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                                             Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), CI.low=round(MDRI.CI[1,1],3), CI.up=round(MDRI.CI[1,2],3)),
                                             Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR[1],3), CI.low=round(FRR.CI[1,1],4), CI.up=round(FRR.CI[1,2],4)),
                                             Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
      } else
        if(BMest=="MDRI.FRR.indep"){ output <- list(Minimum.SS=SS,
                                                    Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)), CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                                                    Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                                                    Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI*365.25,3), CI.low=round(MDRI.CI[,1],3), CI.up=round(MDRI.CI[,2],3)),
                                                    Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR,3), CI.low=round(FRR.CI[,1],4), CI.up=round(FRR.CI[,2],4)),
                                                    Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )

        }



  }

  else
    if(SS=="out" & n2=="out"){
      if(n1=="out"){stop("both n1 and n2 cannot be designated 'out'")}
      if (BMest=="same.test")
      {SS <- ceiling(I[2]^2*((1/PrevH[2])*(DE_H[2]/(1-PrevH[2])+(DE_R[2]/CR[2])*(PrevR[2]/PrevH[2])*(1-PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2]-FRR[2])^2)))/
                       ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)) )^2 -
                          ((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)  -
                          ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2) -
                          1/n1*I[1]^2*((1/PrevH[1])*(DE_H[1]/(1-PrevH[1])+(DE_R[1]/CR[1])*(PrevR[1]/PrevH[1])*(1-PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1]-FRR[1])^2)))
                       )
      )
      } else if(BMest=="FRR.indep")
      { SS <- ceiling(I[2]^2*((1/PrevH[2])*(DE_H[2]/(1-PrevH[2])+(DE_R[2]/CR[2])*(PrevR[2]/PrevH[2])*(1-PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2]-FRR[2])^2))) /
                        ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)))^2  -  ((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2 -
                           sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2) -
                           1/n1*I[1]^2*((1/PrevH[1])*(DE_H[1]/(1-PrevH[1])+(DE_R[1]/CR[1])*(PrevR[1]/PrevH[1])*(1-PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1]-FRR[1])^2)))
                        ))
      } else if(BMest=="MDRI.FRR.indep")
      { SS <- ceiling(I[2]^2*((1/PrevH[2])*(DE_H[2]/(1-PrevH[2])+(DE_R[2]/CR[2])*(PrevR[2]/PrevH[2])*(1-PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2]-FRR[2])^2))) /
                        ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)))^2  -  sum(I^2*(RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2) -
                           sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2 ) -
                           1/n1*I[1]^2*((1/PrevH[1])*(DE_H[1]/(1-PrevH[1])+(DE_R[1]/CR[1])*(PrevR[1]/PrevH[1])*(1-PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1]-FRR[1])^2)))
                        ))
      }

      if(SS<0){stop("No sample size will meet the given constraints")}


      #Now based on derived common SS, output implied summary statistics
      N<-c(n1,SS) #make derived necessary common sample size a vector
      if(sum(PrevR/PrevH*(CR*(N-N*HIV.neg))<10)>0 ) #formula is count of expected recent infections
        warning("Expected count of 'recent' infections is less than 10 for at least one survey")

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
      RSE_deltaI.infSS <-  if(BMest=="MDRI.FRR.indep"){ sqrt(sum(((RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2 +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)*I^2))/deltaI_Est
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
        output <- list(Minimum.SS=SS,
                       Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)),CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                       Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                       Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), CI.low=round(MDRI.CI[1,1],3), CI.up=round(MDRI.CI[1,2],3)),
                       Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR,3), CI.low=round(FRR.CI[,1],4), CI.up=round(FRR.CI[,2],4)),
                       Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
      } else
        if(BMest=="same.test") {output <- list(Minimum.SS=SS,
                                               Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)),CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                                               Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                                               Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), CI.low=round(MDRI.CI[1,1],3), CI.up=round(MDRI.CI[1,2],3)),
                                               Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR[1],3), CI.low=round(FRR.CI[1,1],4), CI.up=round(FRR.CI[1,2],4)),
                                               Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
        } else
          if(BMest=="MDRI.FRR.indep"){ output <- list(Minimum.SS=SS,
                                                      Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)),CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                                                      Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                                                      Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI*365.25,3), CI.low=round(MDRI.CI[,1],3), CI.up=round(MDRI.CI[,2],3)),
                                                      Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR,3), CI.low=round(FRR.CI[,1],4), CI.up=round(FRR.CI[,2],4)),
                                                      Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )

          }



    }

  else
    if(SS=="out" & n1=="both" & n2=="both"){
if (BMest=="same.test")
  {SS <- ceiling(sum(I^2*((1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2))))/
                    ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)) )^2 -
                    ((RSE_MDRI[1]*MDRI[1])/(MDRI[1]-FRR[1]*BigT))^2*(deltaI_Est^2)  -
                    ((FRR[1]*RSE_FRR[1])^2/((MDRI[1]-FRR[1]*BigT)^4)*(PrevH[2]*(MDRI[1]-PrevR[2]/PrevH[2]*BigT)/(1-PrevH[2])-PrevH[1]*(MDRI[1]-PrevR[1]/PrevH[1]*BigT)/(1-PrevH[1]))^2)
                     ))
} else if(BMest=="FRR.indep")
  { SS <- ceiling(sum(I^2*((1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2)))) /
                  ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)))^2  -  ((MDRI[1]*RSE_MDRI[1])^2)*(I[1]/(MDRI[1]-FRR[1]*BigT)-I[2]/(MDRI[1]-FRR[2]*BigT))^2 -
                  sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)
                  ))
} else if(BMest=="MDRI.FRR.indep")
  { SS <- ceiling(sum(I^2*((1/PrevH)*(DE_H/(1-PrevH)+(DE_R/CR)*(PrevR/PrevH)*(1-PrevR/PrevH)/((PrevR/PrevH-FRR)^2)))) /
                      ((deltaI_Est/(qnorm(1-Power)-qnorm(1-alpha/2)))^2  -  sum(I^2*(RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2) -
                      sum(I^2*(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2 )
                      ))
}

  if(SS<0){stop("No sample size will meet the given constraints")}


  #Now based on derived common SS, output implied summary statistics
  N<-c(SS,SS) #make derived necessary common sample size a vector
  if(sum(PrevR/PrevH*(CR*(N-N*HIV.neg))<10)>0 ) #formula is count of expected recent infections
    warning("Expected count of 'recent' infections is less than 10 for at least one survey")

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
  RSE_deltaI.infSS <-  if(BMest=="MDRI.FRR.indep"){ sqrt(sum(((RSE_MDRI*MDRI/(MDRI-FRR*BigT))^2 +(RSE_FRR*FRR*(MDRI-(PrevR/PrevH)*BigT)/((MDRI-FRR*BigT)*(PrevR/PrevH-FRR)))^2)*I^2))/deltaI_Est
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
                   Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)),CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                   Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                   Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), CI.low=round(MDRI.CI[1,1],3), CI.up=round(MDRI.CI[1,2],3)),
                   Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR,3), CI.low=round(FRR.CI[,1],4), CI.up=round(FRR.CI[,2],4)),
                   Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
  } else
    if(BMest=="same.test") {output <- list(Minimum.Common.SS=SS,
                                                Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)),CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                                                Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                                                Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI[1]*365.25,3), CI.low=round(MDRI.CI[1,1],3), CI.up=round(MDRI.CI[1,2],3)),
                                                Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR[1],3), CI.low=round(FRR.CI[1,1],4), CI.up=round(FRR.CI[1,2],4)),
                                                Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
  } else
    if(BMest=="MDRI.FRR.indep"){ output <- list(Minimum.Common.SS=SS,
                                                     Inc.Difference.Statistics=data.frame(deltaI_Est=round(deltaI_Est,3), RSE_deltaI=round(RSE_deltaI,3), RSE_deltaI.infSS=ifelse(RSE_deltaI.infSS<0.001,"<0.001",round(RSE_deltaI.infSS,3)), Power=round(ss.power,3), Power.infSS=ifelse(Power.infSS>0.99,">0.99",round(Power.infSS,3)),CI.low=round(deltaI_CI[1],4),CI.up=round(deltaI_CI[2],4)),
                                                     Implied.Incidence.Statistics=data.frame(Survey=c(1,2), Given.I=round(I,3), RSE_I=round(RSE_I,3), CI.low=round(CI.low,3), CI.up=round(CI.up,3)),
                                                     Implied.MDRI.Statistics=data.frame(Given.MDRI=round(MDRI*365.25,3), CI.low=round(MDRI.CI[,1],3), CI.up=round(MDRI.CI[,2],3)),
                                                     Implied.FRR.Statistics=data.frame(Given.FRR=round(FRR,3), CI.low=round(FRR.CI[,1],4), CI.up=round(FRR.CI[,2],4)),
                                                     Implied.Subject.Counts=round(data.frame(Survey.1=c(HIV.negative=N[1]*HIV.neg[1],HIV.positive=N[1]-N[1]*HIV.neg[1],HIV.post.tested.for.recent=round(CR[1]*(N[1]-N[1]*HIV.neg[1])),Recency.test.pos=round(PrevR[1]/PrevH[1]*(CR[1]*(N[1]-N[1]*HIV.neg[1]) ))),Survey.2=c(HIV.negative=N[2]*HIV.neg[2],HIV.positive=N[2]-N[2]*HIV.neg[2],HIV.post.tested.for.recent=round(CR[2]*(N[2]-N[2]*HIV.neg[2])),Recency.test.pos=round(PrevR[2]/PrevH[2]*(CR[2]*(N[2]-N[2]*HIV.neg[2]) )) )) ) )
  }
   }



  return(output)

  }

## ---- echo=TRUE----------------------------------------------------------
SSPower(I1 = 0.05, I2 = 0.03, PrevH1 = 0.20, PrevH2 = 0.20, n1 = 5000, n2 = 5000, alpha = 0.05, Power = "out", SS = NULL, CR = 1, DE_H = 1, DE_R = 1, BMest = "same.test", MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.20, BigT = 730)

## ---- echo=TRUE----------------------------------------------------------
SSPower(I1 = 0.05, I2 = 0.03, PrevH1 = 0.20, PrevH2 = 0.15, n1 = "both", n2 = "both", alpha = 0.05, Power = 0.8, SS = "out", CR = 1, DE_H = 1, DE_R = 1, BMest = "FRR.indep", MDRI = 200, RSE_MDRI = 0.05, FRR = c(0.01,0.009), RSE_FRR = c(0.20,0.22), BigT = 730)

## ---- echo=TRUE----------------------------------------------------------
SSPower(I1 = 0.05, I2 = 0.03, PrevH1 = 0.20, PrevH2 = 0.15, n1 = 5000, n2 = "out", alpha = 0.05, Power = 0.8, SS = "out", CR = 1, DE_H = 1, DE_R = 1, BMest = "same.test", MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.21, BigT = 730)

## ---- echo=FALSE---------------------------------------------------------
DM_FirstOrderTerms <- function(prevH, prevR, mdri, frr, bigt) {
  fot_prevH <- (prevR - frr)/(((1 - prevH)^2) * (mdri - frr * bigt))  #E.G. d(I)/d(P_H)
  fot_prevR <- prevH/((1 - prevH) * (mdri - frr * bigt))
  fot_mdri <- (frr * prevH - prevR * prevH)/((1 - prevH) * ((mdri - frr * bigt)^2))
  fot_frr <- (prevH * (bigt * prevR - mdri))/((1 - prevH) * ((mdri - frr * bigt)^2))
  return(c(fot_prevH, fot_prevR, fot_mdri, fot_frr))
}

SSCprecision <- function(I, RSE_I, PrevH, CR, MDRI, RSE_MDRI, FRR, RSE_FRR, BigT = 730, DE_H = 1, DE_R = 1, n = "out", step = 5) {
  var_list <- list(I = I, RSE_I = RSE_I, PrevH = PrevH, CR = CR, MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR,
                   BigT = BigT, DE_H = DE_H, DE_R = DE_R, n = n, step = step)

  # CHECK TO MAKE SURE ONLY TWO VARIABLES ARE ALLOWED TO VARY
  max.list <- 0
  for (i in 1:length(var_list)) {
    if (length(var_list[[i]]) > 1) {
      max.list = max.list + 1
    }
    if (max.list > 2) {
      stop("only a maximum of 2 variables are allowed to vary")
    }
  }

  if (sum(var_list == "out") > 1) {
    stop("only one of the variables RSE_I or n can be requested at a time")
  }

  if (length(var_list[1:8]) < 8) {
    stop("Not enough variables have been specified")
  }

  for (i in 1:12) {
    if (length(var_list[[i]]) > 2 | length(var_list[[i]]) < 1) {
      stop(paste("specify (only) min & max values for ", names(var_list)[i]), sep = "")
    }
  }

  for (i in c(1, 3, 4, 6:8)) {
    if (is.numeric(var_list[[1:length(var_list[i])]]) > 0 & (sum(var_list[[i]] <= 1) != length(var_list[[i]]) | sum(var_list[[i]] >=
                                                                                                                    0) != length(var_list[[i]]))) {
      stop("Some input values are less than 0 or greater than 1")
    }
  }

  if (sum(RSE_I != "out") > 0) {
    if (is.numeric(var_list[[1:length(var_list[2])]]) > 0 & (sum(var_list[[2]] <= 1) != length(var_list[[2]]) | sum(var_list[[2]] >=
                                                                                                                    0) != length(var_list[[2]]))) {
      stop("Some input values are less than 0 or greater than 1")
    }
  }

  # above code does what below code does, only in 1 line. Which should we keep?  if (is.numeric(I)>0) {stopifnot (I<=1 &
  # I>=0)} if (is.numeric(RSE_I)>0) {stopifnot (RSE_I<=1 & RSE_I>=0)} if (is.numeric(PrevH)>0) {stopifnot (PrevH<=1 &
  # PrevH>=0)} if (is.numeric(CR)>0) {stopifnot (CR<=1 & CR>=0)} if (is.numeric(MDRI)>0) {stopifnot (MDRI>=0)} if
  # (is.numeric(RSE_MDRI)>0) {stopifnot (RSE_MDRI<=1 & RSE_MDRI>=0)} if (is.numeric(FRR)>0) {stopifnot (FRR<=1 & FRR>=0)}
  # if (is.numeric(RSE_FRR)>0) {stopifnot (RSE_FRR<=1 & RSE_FRR>=0)} if (is.numeric(DE_H)>0) {stopifnot (DE_H>=1)} if
  # (is.numeric(DE_R)>0) {stopifnot (DE_R>=1)} if (is.numeric(n)>0) {stopifnot (n>100)}

  if (sum(BigT <= 182)) {
    warning("BigT is smaller than half a year")
  }
  if (sum(BigT < MDRI) > 0) {
    stop("MDRI cannot be greater than BigT")
  }
  if (sum(RSE_MDRI < 0.01) > 0) {
    warning("RSE of estimated MDRI is less than 1%")
  }
  if (sum(FRR == 0) > 0) {
    warning("Zero estimated FRR")
  }
  if (sum(FRR > 0.1) > 0) {
    warning("Estimated FRR is greater than 10%")
  }
  if (sum(RSE_FRR > 0.3) > 0) {
    warning("RSE of estimated FRR is greater than 30%")
  }
  if (sum(RSE_FRR < 0.05) > 0) {
    warning("RSE of estimated FRR is less than 5%")
  }
  if (sum(I > 0.2) > 0) {
    warning(paste("Possible error in incidence input.", max(I), " seems exceptionally high", sep = ""))
  }
  if (sum(n < 1000) > 0) {
    warning("Sample size is smaller than 1000")
  }


  # Need error that reads: 'ERROR: Test properties not consistent with test for recent infection'


  ######### THIS WHOLE SECTION HERE IS FOR IF ONE OR MORE OF THE VARIABLES IS ALLOWED TO VARY####### CREATES TWO NULL VALUES, I
  ######### BELIVE FOR WHICH VARIABLES ARE TO VARY
  vary1 <- NULL
  vary2 <- NULL
  # NOW FOR EACH VARIABLE IN LIST, IF IT'S ALLOWED TO VARY, THEN MAKE A MATRIX OF VALUES FOR THE STEP, FROM BEGINNING TO
  # END, STEP NUMBER OF COLUMNS. IF THE VARIABLE IS THE FIRST ONE IN THIS LIST TO VARY THEN WE'RE CALLING THE VARIABLE
  # IT'S SAME NAME, AS A MATRIX, AND IF IT'S THE SECOND, WE'RE CALLING THAT VARIABLE AS IT'S SAME NAME, BUT AS A
  # TRANSPOSED MATRIX, SO NOW THE ROWS ARE THE STEPPING VALUES there's got to be a way to make this more efficient, like
  # if two have been met, STOP...
  if (length(I) == 2) {
    I <- matrix(rep(seq(from = min(I), to = max(I), length.out = step), times = step), ncol = step, nrow = step)
    if (length(vary1) == 0) {
      vary1 <- I
      vary_name1 <- "I"
    } else {
      vary2 = I = t(I)
      vary_name2 = "I"
    }
  }
  if (length(RSE_I) == 2) {
    RSE_I <- matrix(rep(seq(from = min(RSE_I), to = max(RSE_I), length.out = step), times = step), ncol = step, nrow = step)
    if (length(vary1) == 0) {
      vary1 <- RSE_I
      vary_name1 <- "I"
    } else {
      vary2 = RSE_I = t(RSE_I)
      vary_name2 = "I"
    }
  }
  if (length(PrevH) == 2) {
    PrevH <- matrix(rep(seq(from = min(PrevH), to = max(PrevH), length.out = step), times = step), ncol = step, nrow = step)
    if (length(vary1) == 0) {
      vary1 <- PrevH
      vary_name1 <- "PrevH"
    } else {
      vary2 = PervH = t(PrevH)
      vary_name2 <- "PrevH"
    }
  }
  if (length(CR) == 2) {
    CR <- matrix(rep(seq(from = min(CR), to = max(CR), length.out = step), times = step), ncol = step, nrow = step)
    if (length(vary1) == 0) {
      vary1 <- CR
      vary_name1 <- "CR"
    } else {
      vary2 = CR = t(CR)
      vary_name2 = "CR"
    }
  }
  if (length(MDRI) == 2) {
    MDRI <- matrix(rep(seq(from = min(MDRI), to = max(MDRI), length.out = step), times = step), ncol = step, nrow = step)
    # IF THIS VARIABLE IS TO VARY, MAKE A MATRIX OF DIM (STEP*STEP) WHERE EACH COLUMN IS THE SEQUENCE OF VALUES
    if (length(vary1) == 0) {
      vary1 <- MDRI
      vary_name1 <- "MDRI"
    } else {
      vary2 = MDRI = t(MDRI)
      vary_name2 = "MDRI"
    }
  }
  if (length(RSE_MDRI) == 2) {
    RSE_MDRI <- matrix(rep(seq(from = min(RSE_MDRI), to = max(RSE_MDRI), length.out = step), times = step), ncol = step,
                       nrow = step)
    # IF THIS VARIABLE IS TO VARY, MAKE A MATRIX OF DIM (STEP*STEP) WHERE EACH COLUMN IS THE SEQUENCE OF VALUES
    if (length(vary1) == 0) {
      vary1 <- RSE_MDRI
      vary_name1 <- "RSE_MDRI"
    } else {
      vary2 = RSE_MDRI = t(RSE_MDRI)
      vary_name2 = "RSE_MDRI"
    }
  }
  if (length(FRR) == 2) {
    FRR <- matrix(rep(seq(from = min(FRR), to = max(FRR), length.out = step), times = step), nrow = step, ncol = step)
    if (length(vary1) == 0) {
      vary1 <- FRR
      vary_name1 <- "FRR"
    } else {
      vary2 = FRR = t(FRR)
      vary_name2 = "FRR"
    }
  }
  if (length(RSE_FRR) == 2) {
    RSE_FRR <- matrix(rep(seq(from = min(RSE_FRR), to = max(RSE_FRR), length.out = step), times = step), nrow = step,
                      ncol = step)
    if (length(vary1) == 0) {
      vary1 <- RSE_FRR
      vary_name1 <- "RSE_FRR"
    } else {
      vary2 = RSE_FRR = t(RSE_FRR)
      vary_name2 = "RSE_FRR"
    }
  }
  if (length(BigT) == 2) {
    BigT <- matrix(rep(seq(from = min(BigT), to = max(BigT), length.out = step), times = step), nrow = step, ncol = step)
    if (length(vary1) == 0) {
      vary1 <- BigT
      vary_name1 <- "BigT"
    } else {
      vary2 = BigT = t(BigT)
      vary_name2 = "BigT"
    }
  }
  if (length(DE_H) == 2) {
    DE_H <- matrix(rep(seq(from = min(DE_H), to = max(DE_H), length.out = step), times = step), nrow = step, ncol = step)
    if (length(vary1) == 0) {
      vary1 <- DE_H
      vary_name1 <- "DE_H"
    } else {
      vary2 = DE_H = t(DE_H)
      vary_name2 = "DE_H"
    }
  }
  if (length(DE_R) == 2) {
    DE_R <- matrix(rep(seq(from = min(DE_R), to = max(DE_R), length.out = step), times = step), nrow = step, ncol = step)
    if (length(vary1) == 0) {
      vary1 <- DE_R
      vary_name1 <- "DE_R"
    } else {
      vary2 = DE_R = t(DE_R)
      vary_name2 = "DE_R"
    }
  }
  if (length(n) == 2) {
    n <- matrix(rep(seq(from = min(n), to = max(n), length.out = step), times = step), nrow = step, ncol = step)
    if (length(vary1) == 0) {
      vary1 <- n
      vary_name1 <- "n"
    } else {
      vary2 = n = t(n)
      vary_name2 = "n"
    }
  }
  ######### END SECTION FOR IF ONE OR MORE OF THE VARIABLES IS ALLOWED TO VARY#######

  # NOW MAKE VARIABLES IN DAYS TO BE IN UNITS YEARS
  if (is.numeric(MDRI)) {
    MDRI <- MDRI/365.25
  }
  if (is.numeric(BigT)) {
    BigT <- BigT/365.25
  }


  PrevR <- ((I * (1 - PrevH) * (MDRI - FRR * BigT))/PrevH + FRR)
  out2 <- PrevHR <- round(PrevH * PrevR, digits = 5)  #Prev.HIV&recent
  out3 <- PrevHnR <- round(PrevH - PrevHR, digits = 5)  #Prev.HIV&nonrecent

  fot <- DM_FirstOrderTerms(PrevH, PrevR, MDRI, FRR, BigT)
  # if the output of each term of DM_FirstOrderTerms is univariate, do one thing, otherwise, do another...
  if (sum(lengths(var_list) > 1) == 0 | length(fot) == 4) {
    fot_PrevH <- fot[1]
    fot_PrevR <- fot[2]
    fot_MDRI <- fot[3]
    fot_FRR <- fot[4]
  } else {
    # this needs to reflect what happens to the output fot matrix which depends on which variables are allowed to vary.
    if (length(I) > 1 & sum(lengths(var_list)) == 14) {
      fot_PrevH <- matrix(fot[1:(step * step)], nrow = step, ncol = step, byrow = FALSE)
      fot_PrevR <- fot[(step * step + 1)]  #things change if and only if only Incidence is allowed to vary.
      fot_MDRI <- matrix(fot[((step * step) + 2):(((step * step) * 2 + 1))], nrow = step, ncol = step)
      fot_FRR <- matrix(fot[(((step * step) * 2 + 2)):length(fot)], nrow = step, ncol = step)
    } else {
      # here is the situation if only a variable besides incidence is allowed to vary, or any two parameters are allowed to
      # vary.
      fot_PrevH <- matrix(fot[1:(step * step)], nrow = step, ncol = step, byrow = FALSE)
      fot_PrevR <- matrix(fot[((step * step) * 1 + 1):(((step * step) * 2))], nrow = step, ncol = step)
      fot_MDRI <- matrix(fot[((step * step) * 2 + 1):(((step * step) * 3))], nrow = step, ncol = step)
      fot_FRR <- matrix(fot[((step * step) * 3 + 1):(((step * step) * 4))], nrow = step, ncol = step)
    }
  }


  # IF SAMPLE SIZE n IS THE OUPUT VARIABLE (SO PRECISION/RSE_I IS FIXED)
  if (sum(n == "out") > 0) {
    out4 <- RSE_I_inf_ss <- round(sqrt((fot_MDRI * RSE_MDRI * MDRI)^2 + (fot_FRR * RSE_FRR * FRR)^2)/I, digits = 5)  #RSE.I.inf.sample
    out1 <- n <- ceiling(((fot_PrevH^2) * PrevH * (1 - PrevH) * DE_H + (fot_PrevR^2) * (PrevR * (1 - PrevR) * DE_R/(CR *
                                                                                                                      PrevH)))/((RSE_I^2 - RSE_I_inf_ss^2) * I^2))

    if (sum(((PrevH * (1 - PrevH))/n) * DE_H <= 0) > 0) {
      stop("no sample size will meet input constraints")
    }

    out5 <- RSE_PrevH <- round(sqrt(((PrevH * (1 - PrevH))/n) * DE_H)/PrevH, digits = 5)  #RSE.PrevH
    out6 <- RSE_PrevR <- round(sqrt(((PrevR * (1 - PrevR))/n * CR * PrevH) * DE_R)/PrevR, digits = 5)  #RSE.PrevR
    out_names <- c("sample.size", "Prev.HIV.and.recent", "Prev.HIV.and.nonrecent", "RSE.I.inf.sample", "RSE.PrevH",
                   "RSE.PrevR")
  }


  # IF precision IS THE OUPUT VARIABLE (SO sample size n IS FIXED)
  if (sum(RSE_I == "out") > 0) {
    out4 <- RSE_I_inf_ss <- round(sqrt((fot_MDRI * RSE_MDRI * MDRI)^2 + (fot_FRR * RSE_FRR * FRR)^2)/I, digits = 5)

    out5 <- RSE_PrevH <- round(sqrt(((PrevH * (1 - PrevH))/n) * DE_H)/PrevH, digits = 5)
    out6 <- RSE_PrevR <- round(sqrt(((PrevR * (1 - PrevR))/n * CR * PrevH) * DE_R)/PrevR, digits = 5)

    out1 <- RSE_I <- round(sqrt(((fot_PrevH^2) * PrevH * (1 - PrevH) * DE_H + (fot_PrevR^2) * (PrevR * (1 - PrevR) *
                                                                                                 DE_R/(CR * PrevH)))/(n * I^2) + RSE_I_inf_ss^2), digits = 5)
    out_names <- c("RSE_I", "Prev.HIV.and.recent", "Prev.HIV.and.nonrecent", "RSE.I.inf.sample", "RSE.PrevH", "RSE.PrevR")
  }
  # if (sum(RSE_I > 0.50) >0) { warning('Implied RSE of incidence is greater than 50%') }


  if (sum(lengths(var_list)) == 15) {
    # if two variables are to vary
    variable.1 <- vector(length = step)
    variable.2 <- vector(length = step)
    for (i in 1:step) {
      variable.1[i] <- paste(vary_name1, "=", vary1[i, 1])
      variable.2[i] <- paste(vary_name2, "=", vary2[1, i])
    }

    if (length(out1) > 1) {
      out1 <- data.frame(out1)
      row.names(out1) <- variable.1
      colnames(out1) <- variable.2
    }
    if (length(out2) > 1) {
      out2 <- data.frame(out2)
      row.names(out2) <- variable.1
      colnames(out2) <- variable.2
    }
    if (length(out3) > 1) {
      out3 <- data.frame(out3)
      row.names(out3) <- variable.1
      colnames(out3) <- variable.2
    }
    if (length(out4) > 1) {
      out4 <- data.frame(out4)
      row.names(out4) <- variable.1
      colnames(out4) <- variable.2
    }
    if (length(out5) > 1) {
      out5 <- data.frame(out5)
      row.names(out5) <- variable.1
      colnames(out5) <- variable.2
    }
    if (length(out6) > 1) {
      out6 <- data.frame(out6)
      row.names(out6) <- variable.1
      colnames(out6) <- variable.2
    }


    # out1<-data.frame(out1) out2<-data.frame(out2) out3<-data.frame(out3) out4<-data.frame(out4) out5<-data.frame(out5)
    # out6<-data.frame(out6) row.names(out1) <-variable.2 colnames(out1) <-variable.1 row.names(out2) <-variable.2
    # colnames(out2) <-variable.1 row.names(out3) <-variable.2 colnames(out3) <-variable.1 row.names(out4) <-variable.2
    # colnames(out4) <-variable.1 row.names(out5) <-variable.2 colnames(out5) <-variable.1 row.names(out6) <-variable.2
    # colnames(out6) <-variable.1
  }


  if (sum(lengths(var_list)) == 14) {
    # if just one variable is allowed to vary
    variable.1 <- vector(length = step)
    for (i in c(1:step)) {
      variable.1[i] <- paste(vary_name1, "=", vary1[i, 1])
    }
    if (length(out1) > 1) {
      out1 <- data.frame(variable.1, out1[, 1])
      names(out1) <- c(vary_name1, out_names[1])
    }
    if (length(out2) > 1) {
      out2 <- data.frame(variable.1, out2[, 1])
      names(out2) <- c(vary_name1, out_names[2])
    }
    if (length(out3) > 1) {
      out3 <- data.frame(variable.1, out3[, 1])
      names(out3) <- c(vary_name1, out_names[3])
    }

    if (length(out4) > 1) {
      out4 <- data.frame(variable.1, out4[, 1])
      names(out4) <- c(vary_name1, out_names[4])
    }
    if (length(out5) > 1) {
      out5 <- data.frame(variable.1, out5[, 1])
      names(out5) <- c(vary_name1, out_names[5])
    }
    if (length(out6) > 1) {
      out6 <- data.frame(variable.1, out6[, 1])
      names(out6) <- c(vary_name1, out_names[6])
    }
  }

  output <- list(out1, out2, out3, out4, out5, out6)
  names(output) <- out_names

  return(output)
}

## ---- echo=TRUE----------------------------------------------------------
SSCprecision(I = 0.015, RSE_I = 0.25, PrevH = 0.2, CR = 1, MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 730, DE_H = 1.1, DE_R = 1, n = 'out')

## ---- echo=TRUE----------------------------------------------------------
SSCprecision(I = c(0.015,0.02), RSE_I = 0.25, PrevH = c(0.10,0.20), CR = 1, MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 700, DE_H = 1, DE_R = 1, n = 'out', step = 3)

## ---- echo=TRUE----------------------------------------------------------
SSCprecision(I = 0.017, RSE_I = 'out', PrevH = c(0.10,0.20), CR = 1, MDRI = 211, RSE_MDRI = 0.05, FRR = 0.009, RSE_FRR = 0.2, BigT = 720, n = 5000, step = 5)

