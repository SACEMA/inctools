#' Initialize Inctools.jl
#'
#' Initializes the Julia environment and loads the Inctools.jl package.
#' This function is called automatically when needed, but can also be called
#' explicitly to control when initialization happens.
#'
#' @param julia_path Optional path to Julia executable. If NULL, uses JuliaCall's default.
#' @export
inctools_setup <- function(julia_path = NULL) {
  if (!is.null(julia_path)) {
    Sys.setenv(JULIA_HOME = julia_path)
  }
  .inctools_julia_setup()
}

#' Estimate Prevalence with Confidence Interval
#'
#' Calculates prevalence and confidence interval from survey data.
#'
#' @param pos Number of positive results
#' @param n Total sample size
#' @param de Design effect (default: 1.0)
#' @param ci Whether to compute confidence interval (default: FALSE)
#' @param alpha Significance level (default: 0.05)
#' @return If ci=FALSE, returns a list with `prev` and `se`. If ci=TRUE, returns
#'   a list with `prev`, `se`, and `CI` (a vector of length 2).
#' @export
#' @examples
#' \dontrun{
#' # Simple prevalence
#' prevalence(100, 1000)
#'
#' # With confidence interval
#' prevalence(100, 1000, ci = TRUE)
#' }
prevalence <- function(pos, n, de = 1.0, ci = FALSE, alpha = 0.05) {
  .inctools_julia_setup()

  # Call Julia function (note: Julia uses α, not alpha)
  result <- JuliaCall::julia_call("prevalence", pos, n, de,
                                   ci = ci, α = alpha)

  return(result)
}

#' Estimate Incidence from Prevalence Proportions
#'
#' Estimates HIV incidence using the Kassanjee method from prevalence data.
#'
#' @param prev Prevalence estimate (proportion)
#' @param sigma_prev Standard error of prevalence
#' @param prevR Recent infection prevalence (proportion)
#' @param sigma_prevR Standard error of recent infection prevalence
#' @param mdri Mean duration of recent infection (in days)
#' @param sigma_mdri Standard error of MDRI (in days)
#' @param frr False recent rate (proportion)
#' @param sigma_frr Standard error of FRR
#' @param covar Covariance between prev and prevR (default: 0.0)
#' @param T Time cutoff (default: 730.5 days = 2 years)
#' @param timeconversion Time conversion factor (default: 365.25 days/year)
#' @param bs Number of bootstrap samples (default: 0, uses Delta method)
#' @param alpha Significance level (default: 0.05)
#' @param per Multiplier for incidence rate (default: 1)
#' @return A list containing:
#'   \item{I}{Incidence estimate}
#'   \item{CI}{Confidence interval (vector of length 2)}
#'   \item{RSE}{Relative standard error}
#' @export
#' @examples
#' \dontrun{
#' # Using Delta method
#' incprops(0.20, 0.015, 0.10, 0.02, 130, 15, 0.01, 0.005)
#'
#' # Using bootstrap
#' incprops(0.20, 0.015, 0.10, 0.02, 130, 15, 0.01, 0.005, bs = 2000)
#'
#' # With covariance
#' incprops(0.20, 0.015, 0.10, 0.02, 130, 15, 0.01, 0.005,
#'          covar = 0.0002, bs = 2000)
#' }
incprops <- function(prev, sigma_prev, prevR, sigma_prevR,
                     mdri, sigma_mdri, frr, sigma_frr,
                     covar = 0.0, T = 730.5, timeconversion = 365.25,
                     bs = 0L, alpha = 0.05, per = 1L) {
  .inctools_julia_setup()

  # Ensure integer types for bs and per
  bs <- as.integer(bs)
  per <- as.integer(per)

  # Call Julia function (note: Julia uses α, not alpha)
  result <- JuliaCall::julia_call("incprops",
                                   prev, sigma_prev, prevR, sigma_prevR,
                                   mdri, sigma_mdri, frr, sigma_frr,
                                   covar = covar, T = T,
                                   timeconversion = timeconversion,
                                   bs = bs, α = alpha, per = per)

  return(result)
}

#' Estimate Incidence from Count Data
#'
#' Estimates HIV incidence using the Kassanjee method from raw count data.
#' This function internally computes prevalences and calls incprops.
#'
#' @param n Total sample size
#' @param npos Number of positive results
#' @param ntestR Number tested for recency
#' @param nR Number of recent infections
#' @param mdri Mean duration of recent infection (in days)
#' @param frr False recent rate (proportion)
#' @param de_npos Design effect for prevalence (default: 1.0)
#' @param de_nR Design effect for recent prevalence (default: 1.0)
#' @param sigma_mdri Standard error of MDRI (default: 0.0)
#' @param sigma_frr Standard error of FRR (default: 0.0)
#' @param covar Covariance between prev and prevR (default: 0.0)
#' @param T Time cutoff (default: 730.5 days = 2 years)
#' @param timeconversion Time conversion factor (default: 365.25 days/year)
#' @param bs Number of bootstrap samples (default: 0, uses Delta method)
#' @param alpha Significance level (default: 0.05)
#' @param per Multiplier for incidence rate (default: 1)
#' @return A list containing:
#'   \item{I}{Incidence estimate}
#'   \item{CI}{Confidence interval (vector of length 2)}
#'   \item{RSE}{Relative standard error}
#' @export
#' @examples
#' \dontrun{
#' # Basic usage
#' inccounts(1000, 200, 180, 20, 130, 0.01)
#'
#' # With bootstrap and uncertainty in MDRI/FRR
#' inccounts(1000, 200, 180, 20, 130, 0.01,
#'           sigma_mdri = 15, sigma_frr = 0.005, bs = 2000)
#'
#' # With covariance
#' inccounts(1000, 200, 180, 20, 130, 0.01,
#'           sigma_mdri = 15, sigma_frr = 0.005,
#'           covar = 0.0002, bs = 2000)
#' }
inccounts <- function(n, npos, ntestR, nR, mdri, frr,
                      de_npos = 1.0, de_nR = 1.0,
                      sigma_mdri = 0.0, sigma_frr = 0.0,
                      covar = 0.0, T = 730.5, timeconversion = 365.25,
                      bs = 0L, alpha = 0.05, per = 1L) {
  .inctools_julia_setup()

  # Ensure integer types
  n <- as.integer(n)
  npos <- as.integer(npos)
  ntestR <- as.integer(ntestR)
  nR <- as.integer(nR)
  bs <- as.integer(bs)
  per <- as.integer(per)

  # Call Julia function (note: Julia uses α, not alpha)
  result <- JuliaCall::julia_call("inccounts",
                                   n, npos, ntestR, nR, mdri, frr,
                                   de_npos = de_npos, de_nR = de_nR,
                                   σ_mdri = sigma_mdri, σ_frr = sigma_frr,
                                   covar = covar, T = T,
                                   timeconversion = timeconversion,
                                   bs = bs, α = alpha, per = per)

  return(result)
}

#' Test Difference in Incidence Between Two Populations
#'
#' Tests whether incidence differs significantly between two populations.
#'
#' @param prev Vector of prevalence estimates for two groups (length 2)
#' @param sigma_prev Vector of standard errors of prevalence (length 2)
#' @param prevR Vector of recent infection prevalence (length 2)
#' @param sigma_prevR Vector of standard errors of recent prevalence (length 2)
#' @param mdri Mean duration of recent infection (in days)
#' @param sigma_mdri Standard error of MDRI (in days)
#' @param frr False recent rate (proportion)
#' @param sigma_frr Standard error of FRR
#' @param covar Vector of covariances for each group (default: c(0.0, 0.0))
#' @param T Time cutoff (default: 730.5 days = 2 years)
#' @param timeconversion Time conversion factor (default: 365.25 days/year)
#' @param bs Number of bootstrap samples (default: 0, uses Delta method)
#' @param alpha Significance level (default: 0.05)
#' @param per Multiplier for incidence rate (default: 1)
#' @return A list containing:
#'   \item{Delta}{Difference in incidence}
#'   \item{CI}{Confidence interval for difference}
#'   \item{p}{p-value for test of difference}
#' @export
#' @examples
#' \dontrun{
#' # Compare two populations
#' incdif(c(0.20, 0.15), c(0.015, 0.012),
#'        c(0.10, 0.08), c(0.02, 0.015),
#'        130, 15, 0.01, 0.005, bs = 2000)
#' }
incdif <- function(prev, sigma_prev, prevR, sigma_prevR,
                   mdri, sigma_mdri, frr, sigma_frr,
                   covar = c(0.0, 0.0), T = 730.5, timeconversion = 365.25,
                   bs = 0L, alpha = 0.05, per = 1L) {
  .inctools_julia_setup()

  # Ensure integer types
  bs <- as.integer(bs)
  per <- as.integer(per)

  # Call Julia function (note: Julia uses Greek letters α, σ)
  result <- JuliaCall::julia_call("incdif",
                                   prev, sigma_prev, prevR, sigma_prevR,
                                   mdri, sigma_mdri, frr, sigma_frr,
                                   covar = covar, T = T,
                                   timeconversion = timeconversion,
                                   bs = bs, α = alpha, per = per)

  return(result)
}
