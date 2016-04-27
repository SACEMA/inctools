# This script was written by Petra Baumler as a test for different ways to make Alex's idea for hybrid
# delta-method/bootstrap confidence intervals faster. We will not use this in the final package.

BS_Count = 10000
# BS_Vars = NULL
BS_Vars = c("PrevH", "PrevR", "MDRI")
N <- c(5000, 5000)
N_H <- c(1000, 1000)
N_testR <- c(1000, 1000)
N_R <- c(100, 70)
MDRI <- 200
RSE_MDRI <- 0.05
FRR <- 0.01
RSE_FRR <- 0.2
BigT <- 730
Covar_HR <- c(0, 0)
BSDM_spread <- 1000
no_s <- length(N)

PrevH <- N_H/N
PrevR <- N_R/N_testR

MDRI <- MDRI/365.25
BigT <- BigT/365.25

I_Est <- function(PrevH, PrevR, MDRI, FRR, BigT) {
  PrevH * (PrevR - FRR)/((1 - PrevH) * (MDRI - FRR * BigT))
}
I_est <- I_Est(PrevH, PrevR, MDRI, FRR, BigT)
deltaI_est <- I_est[1] - I_est[2]

if (is.element("PrevH", BS_Vars)) {
  BS_Var_PrevH <- (PrevH * (1 - PrevH))/N
  DM_Var_PrevH <- c(0, 0)
} else {
  BS_Var_PrevH <- c(0, 0)
  DM_Var_PrevH <- (PrevH * (1 - PrevH))/N
}
if (is.element("PrevR", BS_Vars)) {
  BS_Var_PrevR <- (PrevR * (1 - PrevR))/N_testR
  DM_Var_PrevR <- c(0, 0)
} else {
  BS_Var_PrevR <- c(0, 0)
  DM_Var_PrevR <- (PrevR * (1 - PrevR))/N_testR
}
if (is.element("MDRI", BS_Vars)) {
  BS_Var_MDRI <- (MDRI * RSE_MDRI)^2
  DM_Var_MDRI <- 0
} else {
  BS_Var_MDRI <- 0
  DM_Var_MDRI <- (MDRI * RSE_MDRI)^2
}
if (is.element("FRR", BS_Vars)) {
  BS_Var_FRR <- (FRR * RSE_FRR)^2
  DM_Var_FRR <- 0
} else {
  BS_Var_FRR <- 0
  DM_Var_FRR <- (FRR * RSE_FRR)^2
}

require(MASS)
BS_RootEst <- matrix(nrow = BS_Count, ncol = 8)
for (i in c(1, 2)) {
  Mu <- c(PrevH[i], PrevR[i], MDRI, FRR)
  sigma <- matrix(c(BS_Var_PrevH[i], Covar_HR[i], 0, 0, Covar_HR[i], BS_Var_PrevR[i], 0, 0, 0, 0, BS_Var_MDRI, 0, 0, 0,
                    0, BS_Var_FRR), nrow = 4, ncol = 4)
  BS_RootEst[, (i * 4 - 3):(i * 4)] <- mvrnorm(n = BS_Count, mu = Mu, Sigma = sigma, empirical = TRUE)
}

I_BSVec1 <- I_Est(BS_RootEst[, 1], BS_RootEst[, 2], BS_RootEst[, 3], BS_RootEst[, 4], BigT = BigT)
I_BSVec2 <- I_Est(BS_RootEst[, 5], BS_RootEst[, 6], BS_RootEst[, 7], BS_RootEst[, 8], BigT = BigT)
I_BSVec1_sort <- sort(I_BSVec1, decreasing = FALSE)
I_BSVec2_sort <- sort(I_BSVec2, decreasing = FALSE)
I_BSVecMat <- cbind(I_BSVec1_sort, I_BSVec2_sort)


deltaI_BSVec <- I_BSVecMat[, 1] - I_BSVecMat[, 2]
BS_Var_I <- c(var(I_BSVecMat[, 1]), var(I_BSVecMat[, 2]))
BS_Var_deltaI <- var(deltaI_BSVec)

fot_PrevH <- (PrevR - FRR)/(((1 - PrevH)^2) * (MDRI - FRR * BigT))
fot_PrevR <- PrevH/((1 - PrevH) * (MDRI - FRR * BigT))
fot_MDRI <- (FRR * PrevH - PrevR * PrevH)/((1 - PrevH) * ((MDRI - FRR * BigT)^2))
fot_FRR <- (PrevH * (BigT * PrevR - MDRI))/((1 - PrevH) * ((MDRI - FRR * BigT)^2))

DM_Var_I <- (fot_PrevH^2) * DM_Var_PrevH + (fot_PrevR^2) * DM_Var_PrevR + (fot_MDRI^2) * DM_Var_MDRI + (fot_FRR^2) * DM_Var_FRR
DM_SD_I <- sqrt(DM_Var_I)
DM_Var_deltaI <- sum((fot_PrevH^2) * DM_Var_PrevH) + sum((fot_PrevR^2) * DM_Var_PrevR) + (fot_MDRI[1] - fot_MDRI[2])^2 *
  DM_Var_MDRI + (fot_FRR[1] - fot_FRR[2])^2 * DM_Var_FRR

Var_I <- BS_Var_I + DM_Var_I
RSE_I <- sqrt(Var_I)/I_est
SD_I <- sqrt(Var_I)
Var_deltaI <- DM_Var_deltaI + BS_Var_deltaI
RSE_deltaI <- sqrt(Var_deltaI)/deltaI_est
CIlo_I <- qnorm(0.025, mean = I_est, sd = sqrt(Var_I))
CIup_I <- qnorm(0.975, mean = I_est, sd = sqrt(Var_I))


######################################################################################################
I_BS_spreadbyDM <- function() {
  CIlo_I <- vector(length = no_s)
  CIup_I <- vector(length = no_s)
  mean_Ispread <- vector(length = no_s)
  RSE_Ispread <- vector(length = no_s)
  for (j in c(1:no_s)) {
    i_bsvec <- I_BSVecMat[, j]
    spread_seq <- (rnorm(BSDM_spread, mean(I_BSVecMat[, j]), DM_SD_I[j])) - mean(I_BSVecMat[, j])
    I_BS_spread <- vector(length = BSDM_spread * BS_Count)
    for (i in c(1:BS_Count)) {
      I_BS_spread[(i * BSDM_spread - BSDM_spread + 1):(i * BSDM_spread)] <- spread_seq + i_bsvec[i]
    }
    CIlo_I[j] <- quantile(I_BS_spread, 0.025)
    CIup_I[j] <- quantile(I_BS_spread, 0.975)
    mean_Ispread[j] <- mean(I_BS_spread)
    RSE_Ispread[j] <- sqrt(var(I_BS_spread))/I_est[j]
  }
  output <- data.frame(mean_BSbyDM = mean_Ispread, RSE_BSbyDM = RSE_Ispread, CIlo = CIlo_I, CIup = CIup_I)
  return(output)
}

######################################################################################################## I_spread <- function(I_BSVec, DM_Var_I) { I_BSVec_sort<-sort(I_BSVec, decreasing=FALSE) AUC_i <- function(x){ AUC <-
######################################################################################################## vector(length=BS_Count) for (i in c(1:BS_Count)) { if((I_BSVec_sort[i]>(I_BSVec_sort[x]-4*sqrt(DM_Var_I))) &
######################################################################################################## (I_BSVec_sort[i]<(I_BSVec_sort[x]+4*sqrt(DM_Var_I)))){ AUC[i] <- pnorm(I_BSVec_sort[x], mean=I_BSVec_sort[i],
######################################################################################################## sd=sqrt(DM_Var_I)) - pnorm(I_BSVec_sort[x-1], mean=I_BSVec_sort[i], sd=sqrt(DM_Var_I))} else {AUC[i] <- NA}}
######################################################################################################## return(sum(AUC, na.rm=TRUE))} distvec<-vector(length=BS_Count-1) for (x in c(2:BS_Count)) {distvec[x] <- AUC_i(x)}
######################################################################################################## I_spreadVec <- distvec/sum(distvec) CD <- vector() for (i in c(1:(BS_Count-1))) {CD[i]<-sum(I_spreadVec[c(1:i)])}
######################################################################################################## xi_lo <- findInterval(0.025,CD) xi_up <- findInterval(0.975, CD) CIlo_I <- I_BSVec_sort[xi_lo] CIup_I <-
######################################################################################################## I_BSVec_sort[(xi_up+1)] output<-data.frame('CIlo_Ispread'=CIlo_I, 'CIup_Ispread'=CIup_I) return (output)}

##################### space range of I_est+/-4sigma equally
##################### #####################################################################################
I_spread2 <- function() {
  SpreadMat <- matrix(ncol = no_s, nrow = BSDM_spread)
  for (i in c(1:no_s)) {
    SpreadMat[, i] <- seq(from = (I_est[i] - 5 * sqrt(Var_I[i])), to = (I_est[i] + 5 * sqrt(Var_I[i])), length.out = BSDM_spread)
  }

  AUCMat <- matrix(ncol = no_s, nrow = BSDM_spread)
  for (i in c(1:no_s)) {
    spread_seq <- SpreadMat[, i]
    i_bsvec <- I_BSVecMat[, i]
    dm_sd = DM_SD_I[i]
    auc <- vector(length = BSDM_spread)
    for (j in c(1:BSDM_spread)) {
      auc[j] <- sum(pnorm(spread_seq[j], mean = i_bsvec, sd = dm_sd))
    }
    AUC <- auc/BS_Count
    AUCMat[, i] <- AUC
  }

  CIlo_I <- vector(length = no_s)
  CIup_I <- vector(length = no_s)
  for (i in c(1:no_s)) {
    xi_lo <- findInterval(0.025, AUCMat[, i])
    xi_up <- findInterval(0.975, AUCMat[, i])
    CIlo_I[i] <- SpreadMat[xi_lo, i]
    CIup_I[i] <- SpreadMat[(xi_up + 1), i]
  }
  output <- data.frame(CIlo = CIlo_I, CIup = CIup_I)
  return(output)
}

########################### == space cummulative density of I_est equally
########################### ==#############################################################################
I_spread3 <- function() {

  p_seq <- seq(from = 1e-04, to = 0.9999, length.out = BSDM_spread)

  seq_by_qnorm <- function(i_est, sd_i) {
    spread_seq <- vector(length = BSDM_spread)
    for (i in c(1:BSDM_spread)) {
      spread_seq[i] <- t(qnorm(p_seq[i], mean = i_est, sd = sd_i))
    }
    return(spread_seq)
  }

  SpreadMat <- matrix(ncol = no_s, nrow = BSDM_spread)
  for (i in c(1:no_s)) {
    Spread_Seq <- seq_by_qnorm(i_est = I_est[i], sd_i = SD_I[i])
    SpreadMat[, i] <- Spread_Seq
  }

  AUCMat <- matrix(ncol = no_s, nrow = BSDM_spread)
  for (i in c(1:no_s)) {
    spread_seq <- SpreadMat[, i]
    i_bsvec <- I_BSVecMat[, i]
    var_i <- Var_I[i]
    dm_sd = DM_SD_I[i]
    auc <- vector(length = BSDM_spread)
    for (j in c(1:BSDM_spread)) {
      auc[j] <- sum(pnorm(spread_seq[j], mean = i_bsvec, sd = dm_sd))
    }
    AUC <- auc/BS_Count
    AUCMat[, i] <- AUC
  }

  CIlo_I <- vector(length = no_s)
  CIup_I <- vector(length = no_s)
  for (i in c(1:no_s)) {
    xi_lo <- findInterval(0.025, AUCMat[, 1])
    xi_up <- findInterval(0.975, AUCMat[, 1])
    CIlo_I[i] <- SpreadMat[xi_lo, i]
    CIup_I[i] <- SpreadMat[(xi_up + 1), i]
  }
  output <- data.frame(CIlo = CIlo_I, CIup = CIup_I)
  return(output)
}


##################### spread2 - fast#####################################################################################
I_spread2.2 <- function() {
  SpreadMat <- matrix(ncol = no_s, nrow = (BSDM_spread))
  for (i in c(1:no_s)) {
    SpreadMat[, i] <- seq(from = (I_est[i] - 5 * sqrt(Var_I[i])), to = (I_est[i] + 5 * sqrt(Var_I[i])), length.out = BSDM_spread)
  }
  AUCMat <- matrix(ncol = no_s, nrow = (BSDM_spread - 1))
  for (i in c(1:no_s)) {
    spread_seq <- SpreadMat[, i]
    i_bsvec <- I_BSVecMat[, i]
    dm_sd = DM_SD_I[i]
    auc <- vector(length = (BSDM_spread - 1))
    for (j in c(1:(BSDM_spread - 1))) {
      i_bsvec_use <- i_bsvec[(i_bsvec < (spread_seq[j] + 5 * dm_sd)) & (i_bsvec > (spread_seq[j] - 5 * dm_sd))]
      auc[j] <- sum((pnorm(spread_seq[j + 1], mean = i_bsvec_use, sd = dm_sd)) - (pnorm(spread_seq[j], mean = i_bsvec_use,
                                                                                        sd = dm_sd)))
    }
    AUC <- auc/sum(auc)
    AUCMat[, i] <- AUC
  }

  CIlo_I <- vector(length = no_s)
  CIup_I <- vector(length = no_s)
  for (i in c(1:no_s)) {
    AUC <- AUCMat[, i]
    CD <- vector()
    for (j in c(1:(BSDM_spread - 1))) {
      CD[j] <- sum(AUC[1:j])
    }
    xi_lo <- findInterval(0.025, CD)
    xi_up <- findInterval(0.975, CD)
    CIlo_I[i] <- SpreadMat[xi_lo, i]
    CIup_I[i] <- SpreadMat[(xi_up + 1), i]
  }
  output <- data.frame(CIlo = CIlo_I, CIup = CIup_I)
  return(output)
}

##################### spread3 - fast#####################################################################################

I_spread3.2 <- function() {
  p_seq <- seq(from = 1e-04, to = 0.9999, length.out = BSDM_spread)
  seq_by_qnorm <- function(i_est, sd_i) {
    spread_seq <- vector(length = BSDM_spread)
    for (i in c(1:BSDM_spread)) {
      spread_seq[i] <- t(qnorm(p_seq[i], mean = i_est, sd = sd_i))
    }
    return(spread_seq)
  }

  SpreadMat <- matrix(ncol = no_s, nrow = BSDM_spread)
  for (i in c(1:no_s)) {
    Spread_Seq <- seq_by_qnorm(i_est = I_est[i], sd_i = SD_I[i])
    SpreadMat[, i] <- Spread_Seq
  }

  AUCMat <- matrix(ncol = no_s, nrow = (BSDM_spread - 1))
  for (i in c(1:no_s)) {
    spread_seq <- SpreadMat[, i]
    i_bsvec <- I_BSVecMat[, i]
    dm_sd = DM_SD_I[i]
    auc <- vector(length = (BSDM_spread - 1))
    for (j in c(1:(BSDM_spread - 1))) {
      i_bsvec_use <- i_bsvec[(i_bsvec < (spread_seq[j] + 5 * dm_sd)) & (i_bsvec > (spread_seq[j] - 5 * dm_sd))]
      auc[j] <- sum((pnorm(spread_seq[j + 1], mean = i_bsvec_use, sd = dm_sd)) - (pnorm(spread_seq[j], mean = i_bsvec_use,
                                                                                        sd = dm_sd)))
    }
    AUC <- auc/sum(auc)
    AUCMat[, i] <- AUC
  }

  CIlo_I <- vector(length = no_s)
  CIup_I <- vector(length = no_s)
  for (i in c(1:no_s)) {
    AUC <- AUCMat[, i]
    CD <- vector()
    for (j in c(1:(BSDM_spread - 1))) {
      CD[j] <- sum(AUC[1:j])
    }
    xi_lo <- findInterval(0.025, CD)
    xi_up <- findInterval(0.975, CD)
    CIlo_I[i] <- SpreadMat[xi_lo, i]
    CIup_I[i] <- SpreadMat[(xi_up + 1), i]
  }
  output <- data.frame(CIlo = CIlo_I, CIup = CIup_I)
  return(output)
}

################################################################################################## I_diag <- function(no_int, I_est, Var_I, I_BSVec, DM_Var_I) { no_int <- no_int kk <- no_int-2 delta_i <-
################################################################################################## (8*sd(I_BSVec))/kk ii <- seq(from=(mean(I_BSVec)-4*sd(I_BSVec)-delta_i), to=(mean(I_BSVec)+4*sd(I_BSVec)), by=delta_i)
################################################################################################## CD_I_BS <- ecdf(I_BSVec) AUC_I_BS<- CD_I_BS(ii+delta_i) - CD_I_BS(ii) delta_i <- (8*sqrt(DM_Var_I))/kk ii <-
################################################################################################## seq(from=(-4*sqrt(DM_Var_I)-delta_i), to=(4*sqrt(DM_Var_I)), by=delta_i) AUC_I_DM<-pnorm((ii+delta_i),mean=0,
################################################################################################## sd=sqrt(DM_Var_I))-pnorm((ii),mean=0, sd=sqrt(DM_Var_I)) X <- matrix(rep(AUC_I_BS, times=no_int), ncol=no_int,
################################################################################################## nrow=no_int) Y <- matrix(rep(AUC_I_DM, times=no_int),nrow=no_int, ncol=no_int) AUC_3d <- t(X)*Y M_ur <-
################################################################################################## sum(diag(AUC_3d))/2+sum(diag(AUC_3d[,-1]))/2 M_ll <- sum(diag(AUC_3d))/2+sum(diag(AUC_3d[-1,]))/2 C_ur <-
################################################################################################## sum(diag(AUC_3d[,-(1:kk)]))/2+(AUC_3d[1,no_int])/2 C_ll <- sum(diag(AUC_3d[-(1:kk),]))/2+(AUC_3d[no_int,1])/2
################################################################################################## sumvec_ur <- vector(length=(kk-1)) for (i in c(1:(kk-1))){sumvec_ur[i] <- sum(diag(AUC_3d[,-(1:i)]))/2 +
################################################################################################## sum(diag(AUC_3d[,-(1:(i+1))]))/2} sumvec_ll <- vector(length=(kk-1)) for (i in c(1:(kk-1))){sumvec_ll[i] <-
################################################################################################## sum(diag(AUC_3d[-(1:i),])/2) + sum(diag(AUC_3d[-(1:(i+1)),]))/2} distvec <- c(C_ll, rev(sumvec_ll), M_ll, M_ur,
################################################################################################## sumvec_ur, C_ur) CD <- vector() for (i in c(1:length(distvec))) {CD[i] <- sum(distvec[c(1:i)])} xi_lo <-
################################################################################################## findInterval(0.025,CD) xi_up <- findInterval(0.975, CD) CIlo_I <- I_est -
################################################################################################## 4*sqrt(Var_I)/no_int*(length(distvec)/2-xi_lo) CIup_I <- I_est + 4*sqrt(Var_I)/no_int*(xi_up+1-length(distvec)/2)
################################################################################################## output<-data.frame('CIlo_Idiag'=CIlo_I, 'CIup_Idiag'=CIup_I) return (output)}

##################################################################################################
I_est
RSE_I
CIup_I
CIlo_I


I_BS_spreadbyDM()

# I_spread (I_BSVec=I_BSVecMat[,1], DM_Var_I=DM_Var_I[1])

I_spread2()

I_spread3()

I_spread2.2()

I_spread3.2()



# I_diag (no_int=300, I_est=I_est[1], Var_I=Var_I[1], I_BSVec=I_BSVecMat[,1], DM_Var_I=DM_Var_I[1])
