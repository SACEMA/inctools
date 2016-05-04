# Incidence Estimation Tools
# Copyright (C) 2015-2016, DST/NRF Centre of Excellence in Epidemiological Modelling and Analysis (SACEMA)
# and individual contributors.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Estimate subject-level false-recent rate for a given time cutoff.
#' Each subject with any observations after the time cutoff is assigned a recency status according to the majority
#' of observations for that subject after the cutoff. In the event of exactly half of the observations being
#' classified as recent, the subject contributes a count of 0.5.
#' The function performs an exact binomial test and reports the estimated probability of testing recent after the
#' cutoff, a confidence interval for the proportion, the number of recent results ('successes'),
#' number of subjects ('trials') and the number of data points contributing to the subject-level estimate.
#' @param data A data frame containing variables for subject identifier, time (since detectable infection), and variables with biomarker readings or recency status (to be specified in recency_vars)
#' @param subid_var The variable in the dataframe identifying subjects
#' @param time_var The variable in the dataframe indicating time between 'time zero' (usually detectable infection) and biomarker measurement
#' @param recency_cutoff_time Recency time cut-off ('Big T'). Default=730.5.
#' @param recency_rule Specified rule for defining recent/non-recent outcomes from biomarker data (see Details)
#' @param recency_vars Variables to be used in determining recency outcomes
#' @param recency_params Vector of numeric parameters (e.g. thresholds) for determining recency according to the relevant rule
#' @param alpha Confidence level, default=0.05.
#' @details
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
#' @export
frrcal <- function(data = data, subid_var = "sid", time_var = "time", recency_cutoff_time = 730.5,
                         recency_rule = "binary_data", recency_vars = "recency_status", recency_params = NULL, alpha = 0.05) {
  names(data)[names(data) == subid_var] <- "sid"
  names(data)[names(data) == time_var] <- "time_since_eddi"
  # Keep only observations after the cuttoff
  data <- subset(data, as.numeric(time_since_eddi) > recency_cutoff_time)
  data <- process_data(data = data, subid_var = subid_var, time_var = time_var, recency_vars = recency_vars,
                       inclusion_time_threshold = 1e+06)
  data <- assign_recency_status(data = data, recency_params = recency_params, recency_rule = recency_rule)
  subjectdata <- data.frame(sid = NA, recent = NA)
  for (subjectid in unique(data$sid)) {
    if (sum(data$recency_status[data$sid == subjectid] == 1)/nrow(data[data$sid == subjectid,
                                                                       ]) == 0.5) {
      subjectdata <- rbind(subjectdata, c(subjectid, 0.5))
    }
    if (sum(data$recency_status[data$sid == subjectid] == 1)/nrow(data[data$sid == subjectid,
                                                                       ]) < 0.5) {
      subjectdata <- rbind(subjectdata, c(subjectid, 0))
    }
    if (sum(data$recency_status[data$sid == subjectid] == 1)/nrow(data[data$sid == subjectid,
                                                                       ]) > 0.5) {
      subjectdata <- rbind(subjectdata, c(subjectid, 1))
    }
  }
  subjectdata <- subset(subjectdata, !is.na(sid))
  binomprob <- binom.test(ceiling(sum(subjectdata$recent)), nrow(subjectdata), p = 0, conf.level = 1 -
                            alpha)
  FRR <- data.frame(round(binomprob$estimate[[1]], 4), round(binomprob$conf.int[1], 4), round(binomprob$conf.int[2],
                                                                                              4), alpha, binomprob$statistic, binomprob$parameter[[1]], nrow(data))
  colnames(FRR) <- c("FRRest", "LB", "UB", "alpha", "n_recent", "n_subjects", "n_observations")
  rownames(FRR) <- ""
  return(FRR)
}
