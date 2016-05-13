# Incidence Estimation Tools Copyright (C) 2015-2016, DST/NRF Centre of
# Excellence in Epidemiological Modelling and Analysis (SACEMA) and individual
# contributors.  This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

#' Power and sample size calculation for assay-based incidence estimation
#'
#' @param I1 Predicted incidence of HIV in survey 1.
#' @param I2 Predicted incidence of HIV in survey 2.
#' @param PrevH1 Predicted prevalence of HIV in survey 1.
#' @param PrevH2 Predicted prevalence of HIV in survey 2.
#' @param n1 Sample size for survey 1. If equal sample sizes for both surveys are desired at a given power level, both n1 and n2 must have value 'both', which is the default.
#not-for-release
# If necessary sample size at a given power level for survey 1 is desired and survey 2 has been completed, n1 must be set to 'out' along with SS.
# If necessary sample size at a given power level for survey 2 is desired and survey 1 has been completed, n2 must be set to 'out' along with SS.
#' @param n2 Sample size for survey 2. If equal sample sizes for both surveys are desired at a given power level, both n1 and n2 must have value 'both', which is the default.
#' @param alpha Significance level for test (default alpha=0.05).
#' @param Power Desired power used to calculate a sample size for the surveys. Default is 0.80, meaning the function outputs the necessary sample size to achieve stated power for a test of differences in incidence. If Power is set to 'out', function will return power of detecting a difference in incidences for given sample size inputs.
#' @param CR Coverage rate (0-1).
#' @param SS Sample size. Default is 'out', meaning the function takes a power argument and outputs a common sample size needed to achieve power level for test of differences for incidence. If power is desired for a given sample size, parameter value is irrelevant; however, values for n1 and n2 must be specified.
#' @param DE_H Design effect of HIV prevalence test (vector/integer). If a single value is specified, that value is assumed to be the value for both surveys.
#' @param DE_R Design effect of recency test (vector/integer). If a single value is specified, that value is assumed to be the value for both surveys.
#' @param BMest Biomarker test parameter (MDRI, FRR, and RSE) estimation by one the 3 options 'same.test'(default), 'FRR.indep', 'MDRI.FRR.indep' (string).
#' @param MDRI mean duration of recent infection [days] (vector/integer). If a single value is specified, that value is assumed to be the value for both surveys.
#' @param RSE_MDRI Relative standard error of MDRI [days] (vector/integer). If a single value is specified, that value is assumed to be the value for both surveys.
#' @param FRR False recent rate (vector/integer). If a single value is specified, that value is assumed to be the value for both surveys.
#' @param RSE_FRR Relative standard error of FRR (vector/integer). If a single value is specified, that value is assumed to be the value for both surveys.
#' @param BigT Post-infection time cut-off (days). Default 730. If a single value is specified, that value is assumed to be the value for both surveys.
#' @return Common sample size of two surveys--or the sample size of one survey given the other has already been completed--necessary to achieve a given power level for testing a null hypothesis that the incidence rates are identical between populations; alternatively, the power of said test under a particular sample size scenario. Function also returns implied statistics from input values on paramters, confidence limits, and population counts.
#' @details The package contains long form documentation in the form of vignettes that cover the use of the main fucntions. Use browseVignettes(package="inctools") to access them.
#'
#' @examples
#'incpower(I1 = 0.05, I2 = 0.03, PrevH1 = 0.20, PrevH2 = 0.20,
#'n1 = 5000, n2 = 5000, alpha = 0.05, Power = "out", SS = NULL,
#'DE_H = c(1,1.1), DE_R = 1, BMest = 'same.test', MDRI = 200,
#'RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.20, BigT = 730)
#'
#'
#'incpower(I1 = 0.05, I2 = 0.03, PrevH1 = 0.20, PrevH2 = 0.20,
#'alpha = 0.05, Power = 0.80, SS = "out", DE_H = 1, DE_R = 1,
#'BMest = 'FRR.indep', MDRI = 200, RSE_MDRI = 0.05,
#'FRR = c(0.01,0.009), RSE_FRR = c(0.20,0.21), BigT = 730)
#'
#'
#not-for-release
# #'incpower(I1 = 0.05, I2 = 0.03, PrevH1 = 0.20, PrevH2 = 0.20,
# #'n1 = 5000, n2 = 'out', alpha = 0.05, Power = 0.80, SS ="out",
# #'DE_H = 1, DE_R = 1, BMest = 'MDRI.FRR.indep', MDRI = 200,
# #'RSE_MDRI = c(0.05,0.06), FRR = c(0.01,0.009),
# #'RSE_FRR = c(0.20,0.21), BigT = 730)
#' @export


incpower <- function(I1, I2, PrevH1, PrevH2, n1 = "both", n2 = "both", alpha = 0.05,
    Power = 0.8, SS = "out", CR = 1, DE_H = 1, DE_R = 1, BMest = "same.test", MDRI,
    RSE_MDRI, FRR, RSE_FRR, BigT = 730) {

    ############ Begin warning messages ################
    stopifnot(PrevH1 <= 1 & PrevH1 >= 0)
    stopifnot(PrevH2 <= 1 & PrevH2 >= 0)
    stopifnot(MDRI >= 0)
    stopifnot(RSE_MDRI <= 1 & RSE_MDRI >= 0)
    stopifnot(FRR <= 1 & FRR >= 0)
    stopifnot(RSE_FRR <= 1 & RSE_FRR >= 0)

    if (sum(BMest == c("same.test", "FRR.indep", "MDRI.FRR.indep")) == 0) {
        stop("BMest option must be same.test, FRR.indep, or MDRI.FRR.indep")
    }

    if (length(MDRI) > length(FRR)) {
        stop("number of inputs for MDRI is larger than number of inputs for FRR")
    }

    if (sum(MDRI < 90) > 0) {
        warning("Estimated MDRI less than 90 days")
    }

    if (sum(FRR > 0.1) > 0) {
        warning("Estimated FRR is greater than 10%")
    }
    if (sum(FRR == 0) > 0) {
        warning("Zero estimated FRR")
    }
    if (sum(RSE_FRR > 0.3) > 0) {
        warning("RSE of estimated FRR is greater than 30%")
    }
    if (sum(RSE_FRR < 0.05) > 0) {
        warning("RSE of estimated FRR is less than 5%")
    }
    if (sum(RSE_MDRI < 0.01) > 0) {
        warning("RSE of estimated MDRI is less than 1%")
    }
    if (sum(RSE_MDRI > 0.3) > 0) {
        warning("RSE of estimated MDRI is greater than 30%")
    }

    if (sum(MDRI > BigT) > 0) {
        stop("MDRI cannot be greater than BigT")
    }

    if (sum(BigT <= 182) > 0) {
        warning("BigT is smaller than half a year")
    }

    # dimension of inputs (number of surveys)
    no_s <- 2
    if (length(MDRI) == 1) {
        MDRI <- rep(MDRI, times = no_s)
    } else {
        MDRI = MDRI
    }
    if (length(FRR) == 1) {
        FRR <- rep(FRR, times = no_s)
    } else {
        FRR = FRR
    }
    if (length(RSE_MDRI) == 1) {
        RSE_MDRI <- rep(RSE_MDRI, times = no_s)
    } else {
        RSE_MDRI = RSE_MDRI
    }
    if (length(RSE_FRR) == 1) {
        RSE_FRR <- rep(RSE_FRR, times = no_s)
    } else {
        RSE_FRR = RSE_FRR
    }
    if (length(DE_H) == 1) {
        DE_H <- rep(DE_H, times = no_s)
    } else {
        DE_H = DE_H
    }
    if (length(DE_R) == 1) {
        DE_R <- rep(DE_R, times = no_s)
    } else {
        DE_R = DE_R
    }
    if (length(CR) == 1) {
        CR <- rep(CR, times = no_s)
    } else {
        CR = CR
    }
    if (BMest == "MDRI.FRR.indep") {
        if (length(BigT) == 1) {
            BigT <- rep(BigT, times = no_s)
        } else {
            BigT = BigT
        }
    }

    ############ End warning messages ################




    MDRI <- MDRI/365.25
    BigT <- BigT/365.25
    if (Power == "out")
        N <- c(n1, n2)
    PrevH <- c(PrevH1, PrevH2)
    I <- c(I1, I2)
    deltaI_Est <- I[1] - I[2]

    HIV.neg <- 1 - PrevH
    PrevR <- I * (1 - PrevH) * (MDRI - FRR * BigT) + (FRR * PrevH)

    DM_Var_MDRI <- (MDRI * RSE_MDRI)^2
    DM_Var_FRR <- (FRR * RSE_FRR)^2

    MDRI.CI <- 365.25 * data.frame(CI.low = stats::qnorm(alpha/2, mean = MDRI, sd = sqrt(DM_Var_MDRI)),
        CI.up = stats::qnorm(1 - alpha/2, mean = MDRI, sd = sqrt(DM_Var_MDRI)))
    FRR.CI <- data.frame(CI.low = stats::qnorm(alpha/2, mean = FRR, sd = sqrt(DM_Var_FRR)),
        CI.up = stats::qnorm(1 - alpha/2, mean = FRR, sd = sqrt(DM_Var_FRR)))


    # Instead of trying to get the 'fot' function to work, I explicitly wrote out the
    # commands for the Var[I], Var[D-I]
    if (Power == "out") {
        if (n1 < 1 | n2 < 1)
            stop("Sample size input must be a positive integer")

        Var_I <- I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
            (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)) + (RSE_MDRI * MDRI/(MDRI -
            FRR * BigT))^2 + (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI -
            FRR * BigT) * (PrevR/PrevH - FRR)))^2)
        RSE_I.infss <- sqrt(I^2 * ((RSE_MDRI * MDRI/(MDRI - FRR * BigT))^2 + (RSE_FRR *
            FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH -
            FRR)))^2))/I
        RSE_I <- sqrt(Var_I)/I
        Var_delta_I <- if (BMest == "MDRI.FRR.indep") {
            Var_I[1] + Var_I[2]
        } else if (BMest == "FRR.indep") {
            sum(I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
                (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)))) + ((MDRI[1] * RSE_MDRI[1])^2) *
                (I[1]/(MDRI[1] - FRR[1] * BigT) - I[2]/(MDRI[1] - FRR[2] * BigT))^2 +
                sum(I^2 * (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI -
                  FRR * BigT) * (PrevR/PrevH - FRR)))^2)
        } else if (BMest == "same.test") {
            sum(I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
                (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)))) + ((RSE_MDRI[1] * MDRI[1])/(MDRI[1] -
                FRR[1] * BigT))^2 * (deltaI_Est^2) + ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] -
                FRR[1] * BigT)^4) * (PrevH[2] * (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 -
                PrevH[2]) - PrevH[1] * (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 -
                PrevH[1]))^2)
        }

        RSE_deltaI <- sqrt(Var_delta_I)/abs(deltaI_Est)
        RSE_deltaI.infSS <- if (BMest == "MDRI.FRR.indep") {
            sqrt(sum(((RSE_MDRI * MDRI/(MDRI - FRR * BigT))^2 + (RSE_FRR * FRR *
                (MDRI - (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH -
                FRR)))^2) * I^2))/deltaI_Est
        } else if (BMest == "FRR.indep") {
            sqrt(((MDRI[1] * RSE_MDRI[1])^2) * (I[1]/(MDRI[1] - FRR[1] * BigT) -
                I[2]/(MDRI[1] - FRR[2] * BigT))^2 + sum(I^2 * (RSE_FRR * FRR * (MDRI -
                (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH - FRR)))^2))/deltaI_Est
        } else if (BMest == "same.test") {
            sqrt(((RSE_MDRI[1] * MDRI[1])/(MDRI[1] - FRR[1] * BigT))^2 * (deltaI_Est^2) +
                ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] - FRR[1] * BigT)^4) * (PrevH[2] *
                  (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 - PrevH[2]) - PrevH[1] *
                  (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 - PrevH[1]))^2))/deltaI_Est
        }



        CI.low <- stats::qnorm(alpha/2, mean = I, sd = sqrt(Var_I))
        CI.up <- stats::qnorm(1 - alpha/2, mean = I, sd = sqrt(Var_I))

        deltaI_CI <- NULL
        deltaI_CI[1] <- stats::qnorm(alpha/2, mean = deltaI_Est, sd = sqrt(Var_delta_I))
        deltaI_CI[2] <- stats::qnorm(1 - alpha/2, mean = deltaI_Est, sd = sqrt(Var_delta_I))

        ss.power <- 1 - pnorm(q = stats::qnorm(1 - alpha/2), mean = 1/RSE_deltaI, sd = 1)
        Power.infSS <- 1 - pnorm(q = stats::qnorm(1 - alpha/2), mean = 1/RSE_deltaI.infSS,
            sd = 1)

        # if(ss.power<0.7){warning('Probability of correct inference less than 70%')}
        # this warning is in spreadsheets, but I don't like it. User expected to have an
        # idea about power.


        if (BMest == "FRR.indep") {
            output <- list(Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI[1] *
                365.25, 3), CI.low = round(MDRI.CI[1, 1], 3), CI.up = round(MDRI.CI[1,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR,
                3), CI.low = round(FRR.CI[, 1], 3), CI.up = round(FRR.CI[, 2], 3)),
                Implied.Subject.Counts = data.frame(Survey.1 = c(HIV.negative = N[1] *
                  HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                  (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                  (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                  HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                  (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                  (CR[2] * (N[2] - N[2] * HIV.neg[2]))))))
        } else if (BMest == "same.test") {
            output <- list(Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI[1] *
                365.25, 3), CI.low = round(MDRI.CI[1, 1], 3), CI.up = round(MDRI.CI[1,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR[1],
                3), CI.low = round(FRR.CI[1, 1], 3), CI.up = round(FRR.CI[1, 2],
                3)), Implied.Subject.Counts = data.frame(Survey.1 = c(HIV.negative = N[1] *
                HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                (CR[2] * (N[2] - N[2] * HIV.neg[2]))))))
        } else if (BMest == "MDRI.FRR.indep") {
            output <- list(Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI *
                365.25, 3), CI.low = round(MDRI.CI[, 1], 3), CI.up = round(MDRI.CI[,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR,
                3), CI.low = round(FRR.CI[, 1], 3), CI.up = round(FRR.CI[, 2], 3)),
                Implied.Subject.Counts = data.frame(Survey.1 = c(HIV.negative = round(N[1] *
                  HIV.neg[1]), HIV.positive = round(N[1] - N[1] * HIV.neg[1]), HIV.post.tested.for.recent = round(CR[1] *
                  (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                  (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = round(N[2] *
                  HIV.neg[2]), HIV.positive = round(N[2] - N[2] * HIV.neg[2]), HIV.post.tested.for.recent = round(CR[2] *
                  (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                  (CR[2] * (N[2] - N[2] * HIV.neg[2]))))))
        }
    } else if (SS == "out" & n1 == "out") {
        if (n2 == "out") {
            stop("both n1 and n2 cannot be designated 'out'")
        }
        if (BMest == "same.test") {
            SS <- ceiling(I[1]^2 * ((1/PrevH[1]) * (DE_H[1]/(1 - PrevH[1]) + (DE_R[1]/CR[1]) *
                (PrevR[1]/PrevH[1]) * (1 - PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1] -
                FRR[1])^2)))/((deltaI_Est/(stats::qnorm(1 - Power) - stats::qnorm(1 - alpha/2)))^2 -
                ((RSE_MDRI[1] * MDRI[1])/(MDRI[1] - FRR[1] * BigT))^2 * (deltaI_Est^2) -
                ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] - FRR[1] * BigT)^4) * (PrevH[2] *
                  (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 - PrevH[2]) - PrevH[1] *
                  (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 - PrevH[1]))^2) - 1/n2 *
                I[2]^2 * ((1/PrevH[2]) * (DE_H[2]/(1 - PrevH[2]) + (DE_R[2]/CR[2]) *
                (PrevR[2]/PrevH[2]) * (1 - PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2] -
                FRR[2])^2)))))
        } else if (BMest == "FRR.indep") {
            SS <- ceiling(I[1]^2 * ((1/PrevH[1]) * (DE_H[1]/(1 - PrevH[1]) + (DE_R[1]/CR[1]) *
                (PrevR[1]/PrevH[1]) * (1 - PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1] -
                FRR[1])^2)))/((deltaI_Est/(stats::qnorm(1 - Power) - stats::qnorm(1 - alpha/2)))^2 -
                ((MDRI[1] * RSE_MDRI[1])^2) * (I[1]/(MDRI[1] - FRR[1] * BigT) - I[2]/(MDRI[1] -
                  FRR[2] * BigT))^2 - sum(I^2 * (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) *
                BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH - FRR)))^2) - 1/n2 * I[2]^2 *
                ((1/PrevH[2]) * (DE_H[2]/(1 - PrevH[2]) + (DE_R[2]/CR[2]) * (PrevR[2]/PrevH[2]) *
                  (1 - PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2] - FRR[2])^2)))))
        } else if (BMest == "MDRI.FRR.indep") {
            SS <- ceiling(I[1]^2 * ((1/PrevH[1]) * (DE_H[1]/(1 - PrevH[1]) + (DE_R[1]/CR[1]) *
                (PrevR[1]/PrevH[1]) * (1 - PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1] -
                FRR[1])^2)))/((deltaI_Est/(stats::qnorm(1 - Power) - stats::qnorm(1 - alpha/2)))^2 -
                sum(I^2 * (RSE_MDRI * MDRI/(MDRI - FRR * BigT))^2) - sum(I^2 * (RSE_FRR *
                FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH -
                FRR)))^2) - 1/n2 * I[2]^2 * ((1/PrevH[2]) * (DE_H[2]/(1 - PrevH[2]) +
                (DE_R[2]/CR[2]) * (PrevR[2]/PrevH[2]) * (1 - PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2] -
                  FRR[2])^2)))))
        }

        if (SS < 0) {
            stop("No sample size will meet the given constraints")
        }


        # Now based on derived common SS, output implied summary statistics

        # make derived necessary common sample size a vector
        N <- c(SS, n2)
        # formula is count of expected recent infections
        if (sum(PrevR/PrevH * (CR * (N - N * HIV.neg)) < 10) > 0)
            warning("Expected count of 'recent' infections is less than 10 for at least one survey")

        Var_I <- I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
            (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)) + (RSE_MDRI * MDRI/(MDRI -
            FRR * BigT))^2 + (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI -
            FRR * BigT) * (PrevR/PrevH - FRR)))^2)
        RSE_I <- sqrt(Var_I)/I
        Var_delta_I <- if (BMest == "MDRI.FRR.indep") {
            Var_I[1] + Var_I[2]
        } else if (BMest == "FRR.indep") {
            sum(I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
                (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)))) + ((MDRI[1] * RSE_MDRI[1])^2) *
                (I[1]/(MDRI[1] - FRR[1] * BigT) - I[2]/(MDRI[1] - FRR[2] * BigT))^2 +
                sum(I^2 * (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI -
                  FRR * BigT) * (PrevR/PrevH - FRR)))^2)
        } else if (BMest == "same.test") {
            sum(I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
                (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)))) + ((RSE_MDRI[1] * MDRI[1])/(MDRI[1] -
                FRR[1] * BigT))^2 * (deltaI_Est^2) + ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] -
                FRR[1] * BigT)^4) * (PrevH[2] * (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 -
                PrevH[2]) - PrevH[1] * (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 -
                PrevH[1]))^2)
        }

        RSE_deltaI <- sqrt(Var_delta_I)/abs(deltaI_Est)
        RSE_deltaI.infSS <- if (BMest == "MDRI.FRR.indep") {
            sqrt(sum(((RSE_MDRI * MDRI/(MDRI - FRR * BigT))^2 + (RSE_FRR * FRR *
                (MDRI - (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH -
                FRR)))^2) * I^2))/deltaI_Est
        } else if (BMest == "FRR.indep") {
            sqrt(((MDRI[1] * RSE_MDRI[1])^2) * (I[1]/(MDRI[1] - FRR[1] * BigT) -
                I[2]/(MDRI[1] - FRR[2] * BigT))^2 + sum(I^2 * (RSE_FRR * FRR * (MDRI -
                (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH - FRR)))^2))/deltaI_Est
        } else if (BMest == "same.test") {
            sqrt(((RSE_MDRI[1] * MDRI[1])/(MDRI[1] - FRR[1] * BigT))^2 * (deltaI_Est^2) +
                ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] - FRR[1] * BigT)^4) * (PrevH[2] *
                  (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 - PrevH[2]) - PrevH[1] *
                  (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 - PrevH[1]))^2))/deltaI_Est
        }

        CI.low <- stats::qnorm(alpha/2, mean = I, sd = sqrt(Var_I))
        CI.up <- stats::qnorm(1 - alpha/2, mean = I, sd = sqrt(Var_I))

        deltaI_CI <- NULL
        deltaI_CI[1] <- stats::qnorm(alpha/2, mean = deltaI_Est, sd = sqrt(Var_delta_I))
        deltaI_CI[2] <- stats::qnorm(1 - alpha/2, mean = deltaI_Est, sd = sqrt(Var_delta_I))

        ss.power <- 1 - pnorm(q = stats::qnorm(1 - alpha/2), mean = 1/RSE_deltaI, sd = 1)
        Power.infSS <- 1 - pnorm(q = stats::qnorm(1 - alpha/2), mean = 1/RSE_deltaI.infSS,
            sd = 1)

        if (BMest == "FRR.indep") {
            output <- list(Minimum.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI[1] *
                365.25, 3), CI.low = round(MDRI.CI[1, 1], 3), CI.up = round(MDRI.CI[1,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR,
                3), CI.low = round(FRR.CI[, 1], 4), CI.up = round(FRR.CI[, 2], 4)),
                Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                  HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                  (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                  (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                  HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                  (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                  (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))
        } else if (BMest == "same.test") {
            output <- list(Minimum.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI[1] *
                365.25, 3), CI.low = round(MDRI.CI[1, 1], 3), CI.up = round(MDRI.CI[1,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR[1],
                3), CI.low = round(FRR.CI[1, 1], 4), CI.up = round(FRR.CI[1, 2],
                4)), Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))
        } else if (BMest == "MDRI.FRR.indep") {
            output <- list(Minimum.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI *
                365.25, 3), CI.low = round(MDRI.CI[, 1], 3), CI.up = round(MDRI.CI[,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR,
                3), CI.low = round(FRR.CI[, 1], 4), CI.up = round(FRR.CI[, 2], 4)),
                Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                  HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                  (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                  (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                  HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                  (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                  (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))

        }



    } else if (SS == "out" & n2 == "out") {
        if (n1 == "out") {
            stop("both n1 and n2 cannot be designated 'out'")
        }
        if (BMest == "same.test") {
            SS <- ceiling(I[2]^2 * ((1/PrevH[2]) * (DE_H[2]/(1 - PrevH[2]) + (DE_R[2]/CR[2]) *
                (PrevR[2]/PrevH[2]) * (1 - PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2] -
                FRR[2])^2)))/((deltaI_Est/(stats::qnorm(1 - Power) - stats::qnorm(1 - alpha/2)))^2 -
                ((RSE_MDRI[1] * MDRI[1])/(MDRI[1] - FRR[1] * BigT))^2 * (deltaI_Est^2) -
                ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] - FRR[1] * BigT)^4) * (PrevH[2] *
                  (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 - PrevH[2]) - PrevH[1] *
                  (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 - PrevH[1]))^2) - 1/n1 *
                I[1]^2 * ((1/PrevH[1]) * (DE_H[1]/(1 - PrevH[1]) + (DE_R[1]/CR[1]) *
                (PrevR[1]/PrevH[1]) * (1 - PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1] -
                FRR[1])^2)))))
        } else if (BMest == "FRR.indep") {
            SS <- ceiling(I[2]^2 * ((1/PrevH[2]) * (DE_H[2]/(1 - PrevH[2]) + (DE_R[2]/CR[2]) *
                (PrevR[2]/PrevH[2]) * (1 - PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2] -
                FRR[2])^2)))/((deltaI_Est/(stats::qnorm(1 - Power) - stats::qnorm(1 - alpha/2)))^2 -
                ((MDRI[1] * RSE_MDRI[1])^2) * (I[1]/(MDRI[1] - FRR[1] * BigT) - I[2]/(MDRI[1] -
                  FRR[2] * BigT))^2 - sum(I^2 * (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) *
                BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH - FRR)))^2) - 1/n1 * I[1]^2 *
                ((1/PrevH[1]) * (DE_H[1]/(1 - PrevH[1]) + (DE_R[1]/CR[1]) * (PrevR[1]/PrevH[1]) *
                  (1 - PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1] - FRR[1])^2)))))
        } else if (BMest == "MDRI.FRR.indep") {
            SS <- ceiling(I[2]^2 * ((1/PrevH[2]) * (DE_H[2]/(1 - PrevH[2]) + (DE_R[2]/CR[2]) *
                (PrevR[2]/PrevH[2]) * (1 - PrevR[2]/PrevH[2])/((PrevR[2]/PrevH[2] -
                FRR[2])^2)))/((deltaI_Est/(stats::qnorm(1 - Power) - stats::qnorm(1 - alpha/2)))^2 -
                sum(I^2 * (RSE_MDRI * MDRI/(MDRI - FRR * BigT))^2) - sum(I^2 * (RSE_FRR *
                FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH -
                FRR)))^2) - 1/n1 * I[1]^2 * ((1/PrevH[1]) * (DE_H[1]/(1 - PrevH[1]) +
                (DE_R[1]/CR[1]) * (PrevR[1]/PrevH[1]) * (1 - PrevR[1]/PrevH[1])/((PrevR[1]/PrevH[1] -
                  FRR[1])^2)))))
        }

        if (SS < 0) {
            stop("No sample size will meet the given constraints")
        }


        # Now based on derived common SS, output implied summary statistics

        # make derived necessary common sample size a vector
        N <- c(n1, SS)

        # formula is count of expected recent infections
        if (sum(PrevR/PrevH * (CR * (N - N * HIV.neg)) < 10) > 0) {
            warning("Expected count of 'recent' infections is less than 10 for at least one survey")
        }

        Var_I <- I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
            (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)) + (RSE_MDRI * MDRI/(MDRI -
            FRR * BigT))^2 + (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI -
            FRR * BigT) * (PrevR/PrevH - FRR)))^2)
        RSE_I <- sqrt(Var_I)/I
        Var_delta_I <- if (BMest == "MDRI.FRR.indep") {
            Var_I[1] + Var_I[2]
        } else if (BMest == "FRR.indep") {
            sum(I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
                (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)))) + ((MDRI[1] * RSE_MDRI[1])^2) *
                (I[1]/(MDRI[1] - FRR[1] * BigT) - I[2]/(MDRI[1] - FRR[2] * BigT))^2 +
                sum(I^2 * (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI -
                  FRR * BigT) * (PrevR/PrevH - FRR)))^2)
        } else if (BMest == "same.test") {
            sum(I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
                (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)))) + ((RSE_MDRI[1] * MDRI[1])/(MDRI[1] -
                FRR[1] * BigT))^2 * (deltaI_Est^2) + ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] -
                FRR[1] * BigT)^4) * (PrevH[2] * (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 -
                PrevH[2]) - PrevH[1] * (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 -
                PrevH[1]))^2)
        }

        RSE_deltaI <- sqrt(Var_delta_I)/abs(deltaI_Est)
        RSE_deltaI.infSS <- if (BMest == "MDRI.FRR.indep") {
            sqrt(sum(((RSE_MDRI * MDRI/(MDRI - FRR * BigT))^2 + (RSE_FRR * FRR *
                (MDRI - (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH -
                FRR)))^2) * I^2))/deltaI_Est
        } else if (BMest == "FRR.indep") {
            sqrt(((MDRI[1] * RSE_MDRI[1])^2) * (I[1]/(MDRI[1] - FRR[1] * BigT) -
                I[2]/(MDRI[1] - FRR[2] * BigT))^2 + sum(I^2 * (RSE_FRR * FRR * (MDRI -
                (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH - FRR)))^2))/deltaI_Est
        } else if (BMest == "same.test") {
            sqrt(((RSE_MDRI[1] * MDRI[1])/(MDRI[1] - FRR[1] * BigT))^2 * (deltaI_Est^2) +
                ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] - FRR[1] * BigT)^4) * (PrevH[2] *
                  (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 - PrevH[2]) - PrevH[1] *
                  (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 - PrevH[1]))^2))/deltaI_Est
        }

        CI.low <- stats::qnorm(alpha/2, mean = I, sd = sqrt(Var_I))
        CI.up <- stats::qnorm(1 - alpha/2, mean = I, sd = sqrt(Var_I))

        deltaI_CI <- NULL
        deltaI_CI[1] <- stats::qnorm(alpha/2, mean = deltaI_Est, sd = sqrt(Var_delta_I))
        deltaI_CI[2] <- stats::qnorm(1 - alpha/2, mean = deltaI_Est, sd = sqrt(Var_delta_I))

        ss.power <- 1 - pnorm(q = stats::qnorm(1 - alpha/2), mean = 1/RSE_deltaI, sd = 1)
        Power.infSS <- 1 - pnorm(q = stats::qnorm(1 - alpha/2), mean = 1/RSE_deltaI.infSS,
            sd = 1)

        if (BMest == "FRR.indep") {
            output <- list(Minimum.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI[1] *
                365.25, 3), CI.low = round(MDRI.CI[1, 1], 3), CI.up = round(MDRI.CI[1,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR,
                3), CI.low = round(FRR.CI[, 1], 4), CI.up = round(FRR.CI[, 2], 4)),
                Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                  HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                  (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                  (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                  HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                  (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                  (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))
        } else if (BMest == "same.test") {
            output <- list(Minimum.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI[1] *
                365.25, 3), CI.low = round(MDRI.CI[1, 1], 3), CI.up = round(MDRI.CI[1,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR[1],
                3), CI.low = round(FRR.CI[1, 1], 4), CI.up = round(FRR.CI[1, 2],
                4)), Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))
        } else if (BMest == "MDRI.FRR.indep") {
            output <- list(Minimum.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI *
                365.25, 3), CI.low = round(MDRI.CI[, 1], 3), CI.up = round(MDRI.CI[,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR,
                3), CI.low = round(FRR.CI[, 1], 4), CI.up = round(FRR.CI[, 2], 4)),
                Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                  HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                  (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                  (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                  HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                  (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                  (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))

        }



    } else if (SS == "out" & n1 == "both" & n2 == "both") {
        if (BMest == "same.test") {
            SS <- ceiling(sum(I^2 * ((1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) *
                (PrevR/PrevH) * (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2))))/((deltaI_Est/(stats::qnorm(1 -
                Power) - stats::qnorm(1 - alpha/2)))^2 - ((RSE_MDRI[1] * MDRI[1])/(MDRI[1] -
                FRR[1] * BigT))^2 * (deltaI_Est^2) - ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] -
                FRR[1] * BigT)^4) * (PrevH[2] * (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 -
                PrevH[2]) - PrevH[1] * (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 -
                PrevH[1]))^2)))
        } else if (BMest == "FRR.indep") {
            SS <- ceiling(sum(I^2 * ((1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) *
                (PrevR/PrevH) * (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2))))/((deltaI_Est/(stats::qnorm(1 -
                Power) - stats::qnorm(1 - alpha/2)))^2 - ((MDRI[1] * RSE_MDRI[1])^2) * (I[1]/(MDRI[1] -
                FRR[1] * BigT) - I[2]/(MDRI[1] - FRR[2] * BigT))^2 - sum(I^2 * (RSE_FRR *
                FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH -
                FRR)))^2)))
        } else if (BMest == "MDRI.FRR.indep") {
            SS <- ceiling(sum(I^2 * ((1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) *
                (PrevR/PrevH) * (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2))))/((deltaI_Est/(stats::qnorm(1 -
                Power) - stats::qnorm(1 - alpha/2)))^2 - sum(I^2 * (RSE_MDRI * MDRI/(MDRI -
                FRR * BigT))^2) - sum(I^2 * (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) *
                BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH - FRR)))^2)))
        }

        if (SS < 0) {
            stop("No sample size will meet the given constraints")
        }


        # Now based on derived common SS, output implied summary statistics make derived
        # necessary common sample size a vector
        N <- c(SS, SS)
        # formula is count of expected recent infections
        if (sum(PrevR/PrevH * (CR * (N - N * HIV.neg)) < 10) > 0) {
            warning("Expected count of 'recent' infections is less than 10 for at least one survey")
        }

        Var_I <- I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
            (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)) + (RSE_MDRI * MDRI/(MDRI -
            FRR * BigT))^2 + (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI -
            FRR * BigT) * (PrevR/PrevH - FRR)))^2)
        RSE_I <- sqrt(Var_I)/I
        Var_delta_I <- if (BMest == "MDRI.FRR.indep") {
            Var_I[1] + Var_I[2]
        } else if (BMest == "FRR.indep") {
            sum(I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
                (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)))) + ((MDRI[1] * RSE_MDRI[1])^2) *
                (I[1]/(MDRI[1] - FRR[1] * BigT) - I[2]/(MDRI[1] - FRR[2] * BigT))^2 +
                sum(I^2 * (RSE_FRR * FRR * (MDRI - (PrevR/PrevH) * BigT)/((MDRI -
                  FRR * BigT) * (PrevR/PrevH - FRR)))^2)
        } else if (BMest == "same.test") {
            sum(I^2 * ((1/N) * (1/PrevH) * (DE_H/(1 - PrevH) + (DE_R/CR) * (PrevR/PrevH) *
                (1 - PrevR/PrevH)/((PrevR/PrevH - FRR)^2)))) + ((RSE_MDRI[1] * MDRI[1])/(MDRI[1] -
                FRR[1] * BigT))^2 * (deltaI_Est^2) + ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] -
                FRR[1] * BigT)^4) * (PrevH[2] * (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 -
                PrevH[2]) - PrevH[1] * (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 -
                PrevH[1]))^2)
        }

        RSE_deltaI <- sqrt(Var_delta_I)/abs(deltaI_Est)
        RSE_deltaI.infSS <- if (BMest == "MDRI.FRR.indep") {
            sqrt(sum(((RSE_MDRI * MDRI/(MDRI - FRR * BigT))^2 + (RSE_FRR * FRR *
                (MDRI - (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH -
                FRR)))^2) * I^2))/deltaI_Est
        } else if (BMest == "FRR.indep") {
            sqrt(((MDRI[1] * RSE_MDRI[1])^2) * (I[1]/(MDRI[1] - FRR[1] * BigT) -
                I[2]/(MDRI[1] - FRR[2] * BigT))^2 + sum(I^2 * (RSE_FRR * FRR * (MDRI -
                (PrevR/PrevH) * BigT)/((MDRI - FRR * BigT) * (PrevR/PrevH - FRR)))^2))/deltaI_Est
        } else if (BMest == "same.test") {
            sqrt(((RSE_MDRI[1] * MDRI[1])/(MDRI[1] - FRR[1] * BigT))^2 * (deltaI_Est^2) +
                ((FRR[1] * RSE_FRR[1])^2/((MDRI[1] - FRR[1] * BigT)^4) * (PrevH[2] *
                  (MDRI[1] - PrevR[2]/PrevH[2] * BigT)/(1 - PrevH[2]) - PrevH[1] *
                  (MDRI[1] - PrevR[1]/PrevH[1] * BigT)/(1 - PrevH[1]))^2))/deltaI_Est
        }

        CI.low <- stats::qnorm(alpha/2, mean = I, sd = sqrt(Var_I))
        CI.up <- stats::qnorm(1 - alpha/2, mean = I, sd = sqrt(Var_I))

        deltaI_CI <- NULL
        deltaI_CI[1] <- stats::qnorm(alpha/2, mean = deltaI_Est, sd = sqrt(Var_delta_I))
        deltaI_CI[2] <- stats::qnorm(1 - alpha/2, mean = deltaI_Est, sd = sqrt(Var_delta_I))

        ss.power <- 1 - pnorm(q = stats::qnorm(1 - alpha/2), mean = 1/RSE_deltaI, sd = 1)
        Power.infSS <- 1 - pnorm(q = stats::qnorm(1 - alpha/2), mean = 1/RSE_deltaI.infSS,
            sd = 1)

        if (BMest == "FRR.indep") {
            output <- list(Minimum.Common.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI[1] *
                365.25, 3), CI.low = round(MDRI.CI[1, 1], 3), CI.up = round(MDRI.CI[1,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR,
                3), CI.low = round(FRR.CI[, 1], 4), CI.up = round(FRR.CI[, 2], 4)),
                Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                  HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                  (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                  (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                  HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                  (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                  (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))
        } else if (BMest == "same.test") {
            output <- list(Minimum.Common.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI[1] *
                365.25, 3), CI.low = round(MDRI.CI[1, 1], 3), CI.up = round(MDRI.CI[1,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR[1],
                3), CI.low = round(FRR.CI[1, 1], 4), CI.up = round(FRR.CI[1, 2],
                4)), Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))
        } else if (BMest == "MDRI.FRR.indep") {
            output <- list(Minimum.Common.SS = SS, Inc.Difference.Statistics = data.frame(deltaI_Est = round(deltaI_Est,
                3), RSE_deltaI = round(RSE_deltaI, 3), RSE_deltaI.infSS = ifelse(RSE_deltaI.infSS <
                0.001, "<0.001", round(RSE_deltaI.infSS, 3)), Power = round(ss.power,
                3), Power.infSS = ifelse(Power.infSS > 0.99, ">0.99", round(Power.infSS,
                3)), CI.low = round(deltaI_CI[1], 4), CI.up = round(deltaI_CI[2],
                4)), Implied.Incidence.Statistics = data.frame(Survey = c(1, 2),
                Given.I = round(I, 3), RSE_I = round(RSE_I, 3), CI.low = round(CI.low,
                  3), CI.up = round(CI.up, 3)), Implied.MDRI.Statistics = data.frame(Given.MDRI = round(MDRI *
                365.25, 3), CI.low = round(MDRI.CI[, 1], 3), CI.up = round(MDRI.CI[,
                2], 3)), Implied.FRR.Statistics = data.frame(Given.FRR = round(FRR,
                3), CI.low = round(FRR.CI[, 1], 4), CI.up = round(FRR.CI[, 2], 4)),
                Implied.Subject.Counts = round(data.frame(Survey.1 = c(HIV.negative = N[1] *
                  HIV.neg[1], HIV.positive = N[1] - N[1] * HIV.neg[1], HIV.post.tested.for.recent = round(CR[1] *
                  (N[1] - N[1] * HIV.neg[1])), Recency.test.pos = round(PrevR[1]/PrevH[1] *
                  (CR[1] * (N[1] - N[1] * HIV.neg[1])))), Survey.2 = c(HIV.negative = N[2] *
                  HIV.neg[2], HIV.positive = N[2] - N[2] * HIV.neg[2], HIV.post.tested.for.recent = round(CR[2] *
                  (N[2] - N[2] * HIV.neg[2])), Recency.test.pos = round(PrevR[2]/PrevH[2] *
                  (CR[2] * (N[2] - N[2] * HIV.neg[2])))))))
        }
    }



    return(output)

}
