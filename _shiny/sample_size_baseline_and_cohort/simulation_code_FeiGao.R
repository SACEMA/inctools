#### Data Generation ####
# N: total screening sample size
# para: lambda0, p, q, r, OmegaT, sigmaOmegaT, betaT, sigmabetaT, Time, tau, alpha, beta
# R: risk ratio, default=1
generatedata <- function(N, para, R=1){
	p = para$p; q = para$q; r = para$r; lambda0 = para$lambda0; lambda1 = lambda0*R; tau = para$tau
	OmegaT = para$OmegaT; sigmaOmegaT = para$sigmaOmegaT; betaT = para$betaT; sigmabetaT = para$sigmabetaT; Time = para$Time
	Npos = rbinom(1,N,p); Nneg = N - Npos
	Npos_test = rbinom(1,Npos,q); Nneg_trial = rbinom(1,Nneg,r)
	PR = betaT + lambda0*(1-p)/p*(OmegaT-betaT*Time)
	NR = rbinom(1,Npos_test,PR)
	wbetaT = rnorm(1,betaT,sigmabetaT); wOmegaT = rnorm(1,OmegaT,sigmaOmegaT)
	Nevent = rpois(1,lambda1*tau*Nneg_trial)
	return(list(Npos=Npos, Nneg=Nneg, Npos_test=Npos_test, Nneg_trial=Nneg_trial, NR=NR, Nevent=Nevent, PR=PR, wbetaT = wbetaT, wOmegaT = wOmegaT))
}
#### Setting Parameters ####
# lambda0, p, q, r, OmegaT, sigmaOmegaT, betaT, sigmabetaT, Time, tau, N, alpha, beta
Setting1_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.8){  ##### MOZ #####
	MDRI = 118/365.25; RSE_MDRI = 0.07; FRR = 0.015; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.0101; p = 0.126;
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}
Setting2_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.8){ ##### SA AGYW 14-17 #####
	MDRI = 118/365.25; RSE_MDRI = 0.07; FRR = 0.015; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.047; p = 0.276;
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}
Setting3_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.8){ ##### SA MSM #####
	MDRI = 118/365.25; RSE_MDRI = 0.07; FRR = 0.015; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.125; p = 0.324;
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}
Setting4_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.8){ ##### USA MSM #####
	MDRI = 142/365.25; RSE_MDRI = 0.1; FRR = 0.01; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.0342; p = 0.145; 
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}


Setting5_para <- function(q=1, r=0.8,tau, N=NULL, alpha=0.05, beta=0.8){ ##### New Pop #####
	MDRI = 140/365.25; RSE_MDRI = 0.12; FRR = 0.015; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.063; p = 0.18; 
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}



### lambda0 estimation ###
rec_est <- function(N,Npos,Npos_test,NR,para){
	pest = Npos/N; qest = Npos_test/Npos; Nneg = N - Npos
	betaT = para$betaT; OmegaT = para$OmegaT; Time = para$Time
	lambdaest = (NR/qest - betaT*Npos) / (Nneg*(OmegaT-betaT*Time))
	para_est = list(lambda0=lambdaest, p=pest, q=qest, OmegaT=OmegaT, sigmaOmegaT = para$sigmaOmegaT, betaT=betaT, sigmabetaT = para$sigmabetaT, Time=Time)
	varlog = rec_gamma_R1(para_est,N)$varlog
	Est=lambdaest; SE = sqrt(varlog)*lambdaest;
	CI = c(lambdaest-qnorm(0.975)*SE, lambdaest+qnorm(0.975)*SE)
	CI_log = c(lambdaest * exp(-qnorm(0.975)*sqrt(varlog)), lambdaest * exp(qnorm(0.975)*sqrt(varlog)))
	return(list(Est=Est,SE = SE,CI = CI, CI_log = CI_log))
}

### lambda1 estimation ###
trial_est <- function(N,Nneg,Nneg_trial,Nevent,tau){
	pest = 1-Nneg/N; rest = Nneg_trial/Nneg
	lambdaest = Nevent /Nneg_trial/tau
	para_est = list(p=pest,r=rest,tau=tau)
	varlog = trial_gamma_T1(para_est,lambdaest, N)$varlog
	Est=lambdaest; SE = sqrt(varlog)*lambdaest;
	CI = c(lambdaest - qnorm(0.975)*SE, lambdaest + qnorm(0.975)*SE)
	CI_log = c(lambdaest * exp(-qnorm(0.975)*sqrt(varlog)), lambdaest * exp(qnorm(0.975)*sqrt(varlog)))
	return(list(Est=Est,SE = SE,CI = CI, CI_log = CI_log))
}
### log(R) = loglambda1 - loglambda0 estimation ###
logHR_est<-function(lambda0_est,lambda1_est){
	lambda0 = lambda0_est$Est; lambda0_SE = lambda0_est$SE; loglambda0_SE = lambda0_SE/lambda0
	lambda1 = lambda1_est$Est; lambda1_SE = lambda1_est$SE; loglambda1_SE = lambda1_SE/lambda1
	logratio = log(lambda1) - log(lambda0)
	logratio_SE = sqrt(loglambda0_SE^2 + loglambda1_SE^2)
	return(list(Est=logratio,SE = logratio_SE, CI = c(logratio - qnorm(0.975)*logratio_SE,logratio + qnorm(0.975)*logratio_SE)))
}
### rho = 1-R estimation ###
efficacy_est<-function(lambda0_est,lambda1_est){
	lambda0 = lambda0_est$Est; lambda0_SE = lambda0_est$SE; loglambda0_SE = lambda0_SE/lambda0
	lambda1 = lambda1_est$Est; lambda1_SE = lambda1_est$SE; loglambda1_SE = lambda1_SE/lambda1
	ratio = lambda1/lambda0
	ratio_SE = sqrt(lambda0_SE^2 *lambda1^2/ lambda0^4 + lambda1_SE^2 / lambda0^2)
	logratio_SE = sqrt(loglambda0_SE^2 + loglambda1_SE^2)
	rho = 1 - ratio
	rho_L = 1- (ratio + qnorm(0.975)*ratio_SE)
	rho_U = 1- (ratio - qnorm(0.975)*ratio_SE)
	rho_L_log = 1-ratio*exp(qnorm(0.975)*logratio_SE)
	rho_U_log = 1-ratio*exp(-qnorm(0.975)*logratio_SE)
	return(list(Est=rho,SE = ratio_SE, CI = c(rho_L,rho_U),CI_log = c(rho_L_log,rho_U_log)))
}