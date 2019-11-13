# Incidence Estimation Tools. Copyright (C) 2015-2019, individual contributors
# and Stellenbosch University.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Estimate False-Recent Rate
#'
#' Estimates subject-level false-recent rate (FRR) for a given time cutoff.
#' Each subject with any observations after the time cutoff is assigned a recency status according to the majority
#' of observations for that subject after the cutoff. In the event of exactly half of the observations being
#' classified as recent, the subject contributes a count of 0.5.
#' The function performs an exact binomial test and reports the estimated probability of testing recent after the
#' cutoff, a confidence interval for the proportion, the number of recent results ('successes'),
#' number of subjects ('trials') and the number of data points contributing to the subject-level estimate.
#'
#' @param data A data frame containing variables for subject identifier, time (since detectable infection), and variables with biomarker readings or recency status (to be specified in recency_vars)
#' @param subid_var The variable in the dataframe identifying subjects
#' @param time_var The variable in the dataframe indicating time between 'time zero' (usually detectable infection) and biomarker measurement
#' @param recency_cutoff_time Recency time cut-off ('Big T'). Default=730.5.
#' @param recency_rule Specified rule for defining recent/non-recent outcomes from biomarker data (see Details)
#' @param recency_vars Variables to be used in determining recency outcomes
#' @param recency_params Vector of numeric parameters (e.g. thresholds) for determining recency according to the relevant rule
#' @param alpha Confidence level, default=0.05.
#' @param method Method for computing confidence interval on binomial probability (passed to binom::binom.confint). Default is Clopper-Pearson 'exact' method. Accepted values: `c("exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit")`.
#' @param debug Enable debugging mode (browser)
#' @details The package contains long form documentation in the form of vignettes that cover the use of the main fucntions. Use browseVignettes(package="inctools") to access them.
#'
#' recency_rule: binary_data - supply a binary variable with 1=recent and 0=non-recent in recency_vars.
#'
#' recency_rule: independent_thresholds: supply one threshold variable per biomarker in recency_vars and the relevant
#' thresholds, as well as whether a value below or above each threshold indicates recency in recency_params.
#'
#' recency_params expects a list of pairs of thresholds and thresholdtypes, with zero indicating a reading below
#' the threshold implies recency and 1 that a reading above the threshold implies recency. (Note: two values,
#' a threshold and a thresholdtype per variable must be specified in recency_params. For example, if you specify
#' recency_vars = c('ODn','ViralLoad') you may specify recency_params = c(1.5,0,500,1), meaning that an ODn reading
#' below 1.5 AND a viral load reasing above 500 indicates a recent result. Objects with missing values in its
#'  biomarker readings will be excluded from caculation.
#' @examples
#' frrcal(data=excalibdata,
#'        subid_var = "SubjectID",
#'        time_var = "DaysSinceEDDI",
#'        recency_cutoff_time = 730.5,
#'        recency_rule = "independent_thresholds",
#'        recency_vars = c("Result","VL"),
#'        recency_params = c(10,0,1000,1),
#'        method = "exact",
#'        alpha = 0.05)
#' @export
frrcal <- function(data = NULL,
                   subid_var = NULL,
                   time_var = NULL ,
                   recency_cutoff_time = 730.5,
                   recency_rule = "binary_data",
                   recency_vars = NULL,
                   recency_params = NULL,
                   alpha = 0.05,
                   method = "exact",
                   debug = FALSE) {

  if (debug) {browser()}

  if (is.null(recency_rule)) {
    stop("Please specify a recency rule")
  }

  if (is.null(recency_vars)) {
    stop("Please specify at least one Recency Variable")
  }

  if (recency_rule != "binary_data" & recency_rule != "independent_thresholds") {
    stop("Please specify a valid recency rule")
  }

  if (recency_rule == "binary_data") {
    if (length(recency_vars) > 1) {
      stop("Binary data should be specified in one recency (outcome) variable.")
    }
    # This line was broken. Should we have a similar check?
    # if (!all(data$recency_vars == 0 | data$recency_vars == 1)) {
    #   stop("Input data is not binary")
    # }
  }

  if (recency_rule == "independent_thresholds" & length(recency_vars) != 0.5 *
      length(recency_params)) {
    stop("The number of recency variables must match the number of recency paramaters.")
  }

  if (is.null(subid_var) | is.null(time_var)) {
    stop("Subject identifier and time variables must be specified.")
  }

  if (is.null(data)) {stop("Error: No dataframe provided.")}
  if (is.null(subid_var)) {stop("Error: No subject identifier variable provided.")}
  if (is.null(time_var)) {stop("Error: No time variable provided.")}
  if (is.null(time_var)) {stop("Error: No recency variables provided variable provided.")}

  if (is.null(method)) {stop("Error: Confidence interval method must be specified.")}
  
  if (length(method) != 1) {stop("Error: Exactly one confidence interval method must be specified.")}
  
  if ( !(method %in% c("exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit"))) {
    stop("Confidence interval method must be one of 'exact', 'ac', 'asymptotic', 'wilson', 'prop.test', 'bayes', 'logit', 'cloglog', 'probit'. See help of binom::binom.test() for further details.")
    }

  # check that subject id, time and recency variables exist
  variables <- colnames(data)
  if (sum(variables == subid_var) != 1) {
    stop(paste("There is no column", subid_var, "in the data frame."))
  }
  if (sum(variables == time_var) != 1) {
    stop(paste("There is no column", time_var, "in the data frame."))
  }
  for (i in 1:length(recency_vars)) {
    if (sum(variables == recency_vars[i]) != 1) {
      stop(paste("There is no column", recency_vars[i], "in the data frame."))
    }
  }

  data <- data %>%
    process_data(subid_var = subid_var,
                 time_var = time_var,
                 recency_vars = recency_vars,
                 inclusion_time_threshold = 1e+06,
                 debug = debug) %>%
    dplyr::filter(.data$time_since_eddi > recency_cutoff_time) %>%
    assign_recency_status(recency_params = recency_params,
                          recency_rule = recency_rule,
                          debug = debug)

  subjectdata <- tibble::tibble(sid = NA, recent = NA)
  for (subjectid in unique(data$sid)) {
    if (sum(data$recency_status[data$sid == subjectid] == 1)/nrow(data[data$sid ==
                                                                       subjectid, ]) == 0.5) {
      subjectdata <- rbind(subjectdata, c(subjectid, 0.5))
    }
    if (sum(data$recency_status[data$sid == subjectid] == 1)/nrow(data[data$sid ==
                                                                       subjectid, ]) < 0.5) {
      subjectdata <- rbind(subjectdata, c(subjectid, 0))
    }
    if (sum(data$recency_status[data$sid == subjectid] == 1)/nrow(data[data$sid ==
                                                                       subjectid, ]) > 0.5) {
      subjectdata <- rbind(subjectdata, c(subjectid, 1))
    }
  }
  subjectdata <- subjectdata[!is.na(subjectdata$sid),]
  nr <- as.integer(ceiling(sum(subjectdata$recent)))
  n <- nrow(subjectdata)
  p <- nr / n
  sigma <- sqrt( (p * (1 - p)) / n )
  binom_ci <- binom::binom.confint(nr, n, conf.level = 1 - alpha, methods = method)
  FRR <- tibble::tibble(FRRest = p, SE = sigma, CI_LB = binom_ci$lower, CI_UB = binom_ci$upper, 
                        alpha = alpha, n_recent = nr, n_subjects = n, n_observations = nrow(data), 
                        ci_method = method)
  options(pillar.sigfig = 6)
  return(FRR)
}
