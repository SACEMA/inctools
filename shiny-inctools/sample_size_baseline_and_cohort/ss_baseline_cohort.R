power <- function(X1, sigma1, X2, sigma2, alpha = 0.05) {
  DeltaTrue <- X2-X1
  SigmaDelta <- sqrt(sigma1^2 + sigma2^2)
  DeltaCrit <- qnorm(p = alpha/2, mean = 0, sd = SigmaDelta, lower.tail = DeltaTrue<0)
  power <- pnorm(DeltaCrit, mean = DeltaTrue, sd = SigmaDelta, lower.tail = DeltaTrue<0)
  return(power)
}

sigma1 <- function(n, Inc, Prev, CR, MDRI, RSE_MDRI, FRR, RSE_FRR, BigT, DE_H, DE_R) {
  RSE <- suppressWarnings(inctools::incprecision(RSE_I = "out", I = Inc, PrevH = Prev, CR = CR, MDRI = MDRI, 
                                                 RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT,
                                                 DE_H = DE_H, DE_R = DE_R, n = n))$RSE_I[[1]]
  sigma1 <- RSE*Inc
  return(sigma1)
}

sigma2 <- function(n, Inc, Prev, FracIncRed, FUT, CohortCR, DE_C) {
  binomialP <- (1-exp(-Inc*(1-FracIncRed)*FUT))
  n_cohort <- (1-Prev)*CohortCR*n
  sigma_SCs <- sqrt(binomialP*(1-binomialP)*n_cohort)
  RSE <- sigma_SCs/(binomialP*n_cohort)
  sigma2 <- RSE*Inc*(1-FracIncRed)*sqrt(DE_C)
  return(sigma2)
}

powerdif <- function(n, Inc, Prev, CR, MDRI, RSE_MDRI, FRR, RSE_FRR, DE_H, DE_R, BigT, FracIncRed, FUT, CohortCR, Power, alpha, DE_C) {
  sigma1 <- sigma1(n=n, Inc = Inc, Prev = Prev, CR = CR, MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT,
                   DE_H = DE_H, DE_R = DE_R)
  sigma2 <- sigma2(n=n, Inc=Inc, Prev = Prev, FracIncRed=FracIncRed, FUT=FUT, CohortCR=CohortCR, DE_C = DE_C)
  PowerN <- power(X1 = Inc, sigma1 = sigma1, X2 = Inc*(1-FracIncRed), sigma2 = sigma2, alpha = alpha) 
  return(Power - PowerN)
}

ss_baseline_cohort <- function(Inc, Prev, FracIncRed, Power = 0.8, alpha = 0.05, MDRI, RSE_MDRI, FRR, RSE_FRR, CR = 1, DE_H = 1, DE_R = 1, DE_C = 1, 
                       BigT = 730, CohortCR = 1, FUT = 1) {
  #browser()

  
    Ns <- seq(1,50000,100)
    powers <- vector(length = length(Ns))
    for (i in 1:length(Ns)) {
      N_tmp <- Ns[i]
      sigma1 <- sigma1(n=N_tmp, Inc = Inc, Prev = Prev, CR = CR, MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT,
                       DE_H = DE_H, DE_R = DE_R)
      sigma2 <- sigma2(n = N_tmp, Inc = Inc, Prev = Prev, FracIncRed = FracIncRed, FUT = FUT, CohortCR = CohortCR, DE_C = DE_C)
      powers[i] <- power(X1 = Inc, sigma1 = sigma1, X2 = Inc*(1-FracIncRed), sigma2 = sigma2, alpha = alpha) 
    }
    plot <- ggplot2::qplot(Ns, powers, geom = "line", xlab = "N", ylab = "Power") +
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1))
    N <- tryCatch({
      stats::uniroot(powerdif, Inc = Inc, Prev = Prev, CR = CR, MDRI = MDRI, RSE_MDRI = RSE_MDRI, 
                     FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT, DE_H = DE_H, DE_R = DE_R, 
                     FracIncRed = FracIncRed, FUT = FUT, CohortCR = CohortCR, alpha = alpha, 
                     Power = Power, DE_C = DE_C, interval = c(1,1000000))$root
    },
    error = function(err) {
      print("Desired power level cannot be achieved!")
      return(NA)
    })
    return(list(RequiredN = round(N), plot = plot))
}