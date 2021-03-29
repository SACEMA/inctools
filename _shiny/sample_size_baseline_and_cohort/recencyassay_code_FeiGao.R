#### Sample Size Calculation based on Z ####
# gamma0 = rec_gamma_R1(para); gamma1 = trial_gamma_T1(para,lambda1)
samplesize_L_varZ <- function(gamma0, gamma1, para, R, R0 = 1){
	alpha = para$alpha; beta = para$beta
	gamma00 = gamma0$gamma00; gamma01 = gamma0$gamma01; gamma1 = gamma1$gamma1
	vZ = varZ(para, R, R0)
	Zsum = qnorm(1-alpha/2) + qnorm(beta,0,sqrt(vZ))
	return((gamma00+gamma1)/(((log(R)-log(R0))/Zsum)^2-gamma01))
}

#### Sample Size Calculation based on Z tilde ####
# gamma0 = rec_gamma_R2(para); gamma1 = trial_gamma_T2(para,lambda1)
samplesize_I <- function(gamma0, gamma1, para, R, R0 = 1){
	alpha = para$alpha; beta = para$beta; lambda0 = para$lambda0
	gamma00 = gamma0$gamma00; gamma01 = gamma0$gamma01; gamma1 = gamma1$gamma1
	Zsum = qnorm(1-alpha/2) + qnorm(beta)
	return((R0^2*gamma00+R^2*gamma1)/(((R-R0)/Zsum)^2-R0^2*gamma01))
}

#### Rejection Indicator based on Z ####
Rejection <- function(lambda0_est, lambda1_est, R0 = 1, Zcut = qnorm(0.975)){
	if (lambda0_est$Est<=0) return(F)
	if (lambda1_est$Est==0) {
		Zstat = (log(lambda0_est$Est))/(lambda0_est$SE/lambda0_est$Est)
		return(abs(Zstat)>Zcut)
	}
	Zstat = efficacy_test(lambda0_est,lambda1_est, R0)
	return(abs(Zstat)>Zcut)
}

#### Rejection Indicator based on Z tilde ####
Rejection_I <- function(lambda0_est, lambda1_est, R0 = 1, Zcut = qnorm(0.975)){
	if (lambda0_est$Est<=0) return(F)
	if (lambda1_est$Est==0) {
		Zstat = lambda0_est$Est/lambda0_est$SE
		return(abs(Zstat)>Zcut)
	}
	Zstat = efficacy_test_I(lambda0_est, lambda1_est, R0)
	return(abs(Zstat)>Zcut)
}

### Z calculation ###
efficacy_test<-function(lambda0_est,lambda1_est, R0 = 1){
	lambda0 = lambda0_est$Est; lambda0_SE = lambda0_est$SE; V0 = (lambda0_SE/lambda0)^2
	lambda1 = lambda1_est$Est; lambda1_SE = lambda1_est$SE; V1 = (lambda1_SE/lambda1)^2
	wR = lambda1/lambda0; V = V0+V1
	Z = (log(wR)-log(R0))/sqrt(V)
	return(Z)
}

### Z tilde calculation ###
efficacy_test_I<-function(lambda0_est,lambda1_est, R0 = 1){
	wlambda0 = lambda0_est$Est; lambda0_SE = lambda0_est$SE; V0 = (lambda0_SE/wlambda0)^2
	lambda1 = lambda1_est$Est; lambda1_SE = lambda1_est$SE; V1 = (lambda1_SE/lambda1)^2
	wR = lambda1/wlambda0; V = R0^2*wlambda0^2*V0+lambda1^2*V1
	Z = (lambda1 - R0* wlambda0)/sqrt(V)
	return(Z)
}

#### Expected Numbers - Recency Assay ####
rec_num <- function(N,para){
	lambda0 = para$lambda0; p = para$p; q = para$q
	OmegaT = para$OmegaT; betaT = para$betaT; Time = para$Time
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	ENR = N*p*PR*q
	return(ENR)
}
#### Expected Numbers - Active Arm Trial ####
inc_num <- function(N,lambda1,para){
	p = para$p; r = para$r; tau = para$tau
	Event = N * tau * (1-p) * r * lambda1
	return(Event)	
}

#### Expected Numbers - Summary ####
summarize_expectnum <-function(N, para, R){
	p = para$p; q = para$q; lambda0 = para$lambda0; MDRI = para$OmegaT; FRR = para$betaT; Time = para$Time;
	r = para$r; tau = para$tau
	return(c(N,N*p*q,rec_num(N,para),N*(1-p)*r*tau,inc_num(N,lambda0*R,para)))	
}

rec_gamma_R1<-function(para, N=NULL){
	lambda0 = para$lambda0; p = para$p; q = para$q; OmegaT = para$OmegaT; betaT = para$betaT; Time = para$Time
	sigma2Omega = para$sigmaOmegaT^2;sigma2beta = para$sigmabetaT^2
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	gamma00 = (PR*(1-PR)/(PR-betaT)^2/q + 1/(1-p) + (1-p*q)*sigma2beta/(PR-betaT)^2/q)/p
	gamma01 = sigma2Omega/(OmegaT-betaT*Time)^2 + sigma2beta*((OmegaT-PR*Time)/(PR-betaT)/(OmegaT-betaT*Time))^2
	return(list(gamma00 = gamma00, gamma01 = gamma01, varlog = gamma00/N+gamma01))
}

rec_gamma_R2<-function(para, N=NULL){
	lambda0 = para$lambda0; p = para$p; q = para$q; OmegaT = para$OmegaT; betaT = para$betaT; Time = para$Time
	sigma2Omega = para$sigmaOmegaT^2;sigma2beta = para$sigmabetaT^2
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	gamma00 = (PR*(1-PR)/(PR-betaT)^2/q + 1/(1-p))/p # neglect last term
	gamma01 = sigma2Omega/(OmegaT-betaT*Time)^2 + sigma2beta*((OmegaT-PR*Time)/(PR-betaT)/(OmegaT-betaT*Time))^2
	return(list(gamma00 = gamma00, gamma01 = gamma01, varlog = gamma00/N+gamma01))
}

trial_gamma_T1<-function(para,lambda1, N = NULL){
	p = para$p; r = para$r; tau = para$tau
	gamma1 = 1 / (lambda1*(1-p)*r*tau)
	return(list(gamma1 = gamma1, varlog = gamma1/N))
}

trial_gamma_T2<-function(para,lambda1, N = NULL){
	p = para$p; r = para$r; tau = para$tau
	gamma1 = exp(-lambda1*tau) / ((1-p)*r*(1-exp(-lambda1*tau)))
	return(list(gamma1 = gamma1, varlog = gamma1/N))
}

#### NR - N+test*wbetaT, N+test, N+, wOmega -wbeta*T, Nevent, N-enroll, NR
varW <-function(p,q,PR,OmegaT,betaT,r,lambda1,tau){
	variance = matrix(0,7,7)
	variance[1,1] = p*q*PR*(1-PR) + p*q*(1-p*q)*(PR-betaT)^2
	variance[2,2] = p*q*(1-p*q)
	variance[3,3] = p*(1-p)
	variance[5,5] = (1-p)*r*lambda1*tau*(1+lambda1*tau*(1-r)+lambda1*tau*r*p)
	variance[6,6] = (1-p)*r*(1-r+p*r)
	variance[7,7] = p*q*PR*(1-PR*q*p)
	variance[1,2] = variance[2,1] = p*q*(1-p*q)*(PR-betaT)
	variance[1,3] = variance[3,1] = p*(1-p)*q*(PR-betaT)
	variance[1,5] = variance[5,1] = -p*(1-p)*(PR-betaT)*q*r*lambda1*tau
	variance[1,6] = variance[6,1] = -p*(1-p)*(PR-betaT)*q*r
	variance[1,7] = variance[7,1] = p*q*PR*(1-PR) + p*q*(1-p*q)*(PR-betaT)*PR
	variance[2,3] = variance[3,2] = p*(1-p)*q
	variance[2,5] = variance[5,2] = -p*(1-p)*q*r*lambda1*tau
	variance[2,6] = variance[6,2] = -p*(1-p)*q*r
	variance[2,7] = variance[7,2] = p*q*(1-p*q)*PR
	variance[3,5] = variance[5,3] = -p*(1-p)*r*lambda1*tau
	variance[3,6] = variance[6,3] = -p*(1-p)*r
	variance[3,7] = variance[7,3] = p*(1-p)*q*PR
	variance[5,6] = variance[6,5] = (1-p)*r*(1-r+p*r)*lambda1*tau
	variance[5,7] = variance[7,5] = -p*(1-p)*PR*q*r*lambda1*tau
	variance[6,7] = variance[7,6] = -p*(1-p)*PR*q*r
	return(variance)
}
# p=para$p; q= para$q; OmegaT = para$OmegaT;betaT=0; sigma2omega = 0; sigma2beta=0; r=para$r;tau =para$tau
var_AB <- function(p,q,PR,OmegaT, betaT,r,lambda1,tau){
	a = matrix(0,7,2)
	a[1,1] = -1/(p*q*(PR-betaT))
	a[2,1] = 1/(p*q)
	a[3,1] = -1/(p*(1-p))
	a[4,1] = 1/OmegaT
	a[5,1] = 1/((1-p)*r*lambda1*tau)
	a[6,1] = -1/((1-p)*r)
	a[1,2] = -2*PR*(1-PR)/(p*q)^2/(PR-betaT)^3
	a[2,2] = PR^2/(p*q*(PR-betaT))^2
	a[3,2] = -(p^2+(1-p)^2)/(p*(1-p))^2
	a[4,2] = 0
	a[5,2] = -1/((1-p)*r*lambda1*tau)^2
	a[6,2] = 0
	a[7,2] = (1-2*PR)/(p*q*(PR-betaT))^2
	cov_mat = varW(p,q,PR,OmegaT, betaT, r,lambda1,tau)
	return(t(a)%*%cov_mat%*%a)
}

varZ <-function(para, R, R0 = 1){
	p=para$p; q= para$q; OmegaT = para$OmegaT;betaT=para$betaT; r=para$r;tau =para$tau; Time = para$Time
	lambda0 = para$lambda0; lambda1 = lambda0*R
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	varAB = var_AB(p,q,PR,OmegaT,betaT,r,lambda1,tau)
	tA = log(lambda1/lambda0)-log(R0)
	tB = PR*(1-PR)/(p*q*(PR-betaT)^2) + 1/p + 1/(1-p) + 1/((1-p)*r*lambda1*tau)
	c = rep(0,2); c[1] = 1/sqrt(tB); c[2] = -tA/(2*tB^(3/2))
	return(t(c)%*%varAB%*%c)
}