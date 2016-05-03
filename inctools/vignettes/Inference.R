## ---- echo=FALSE---------------------------------------------------------
prevcounts <- function(N, N_H, N_testR, N_R, DE_H=1, DE_R=1) {
  if (sum(N_H>N)>0 | sum(N_R>N)>0) {
    stop("sample subset larger than total sample")
  }
  if (sum(N_testR<N_R)>0) {
    stop("counts of recency tested less than counts of those found to be recently infected")
  }
  if (sum(DE_H<1) >0 | sum(DE_R<1) >0) {
    stop("Design effects must be >=1")
  }
  if(sum(N_testR/N_H<.95)>0){
    warning("Less than 95% of HIV-positive subjects have recency test results available")
  }

  no_s <- length(N)
  if (length(DE_H)==1)     {DE_H <- rep(DE_H, times=no_s)}         else {DE_H=DE_H}
  if (length(DE_R)==1)     {DE_R <- rep(DE_R, times=no_s)}         else {DE_R=DE_R}
  stopifnot(no_s==length(N_H)  & no_s==length(N_R) & no_s==length(N_testR) &
            no_s==length(DE_H) & no_s==length(DE_R))
  PrevH <- N_H/N
  PrevR <- N_R/N_testR
  RSE_PrevH <- sqrt(((PrevH*(1-PrevH))/N)*DE_H)/PrevH
  RSE_PrevR <- sqrt(((PrevR*(1-PrevR))/N_testR)*DE_R)/PrevR

  output <- data.frame("PrevH"=PrevH, "PrevR"=PrevR, "RSE_PrevH"=RSE_PrevH, "RSE_PrevR"=RSE_PrevR)
  return(output)
}

## ---- fig.show='hold'----------------------------------------------------
prevcounts(N = 5000, N_H = 1000, N_testR = 1000, N_R = 70, DE_R = 1.1)

## ---- fig.show='hold'----------------------------------------------------
prevcounts (N = c(5000,5000), N_H = c(1000,1000), N_testR = c(1000,1000),
N_R = c(100,70), DE_H = 1, DE_R = c(1,1.1))

