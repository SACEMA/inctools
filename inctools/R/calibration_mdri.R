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

#'Estimate MDRI (point estimate and confidence interval) using binomial regression and a maximum likelihood approach
#'
#' @param data A data frame containing variables for subject identifier, time (since detectable infection), and variables with biomarker readings or recency status (to be specified in recency_vars)
#' @param functional_forms Select functional form/link function combinations for fitting probability of testing recent as a function of time to data using binomial regression
#' (see Details). Default=all supported functional forms.
#' @param subid_var The variable in the dataframe identifying subjects
#' @param time_var The variable in the dataframe indicating time between 'time zero' (usually detectable infection) and biomarker measurement
# @param infwind_var The (optional) variable indicating the size of the
# seroconversion window, in the event time_var is relative to first positive test
# preceded by an infection window.  WE NEED A WINDOW PARAMATER AND WINDOW OPTIONS
# - e.g.  (1) uniform window width or (2) normal assumption, SDs
#' @param recency_rule Specified rule for defining recent/non-recent outcomes from biomarker data (see Details)
#' @param recency_vars Variables to be used in determining recency outcomes
#' @param recency_params Vector of numeric parameters (e.g. thresholds) for determining recency according to the relevant rule
#' @param recency_cutoff_time Recency time cut-off ('Big T'). Default=730.5.
#' @param inclusion_time_threshold Data points beyond this time are excluded from the calculation (in same unit as recency_cutoff_time, default=800).
#' @param n_bootstraps Number of subject-level bootstrap resampling operations for estimating confidence intervals, default=100 (useful for testing purposes only)
#' @param alpha Confidence level, default=0.05.
# ADD OPTION TO GET FULL LIST OF MDRIs from the bootstrapping procedure or the
# shape of the distribution or something
#' @param plot Specifies whether a plot of the probability of testing recent over time should be produced
#' @param parallel Set to TRUE in order to perform bootstrapping in parallel on a multicore or multiprocessor syste. Not available on Windows.
#' @param cores Set number of cores for parallel processing when parallel=TRUE. This defaults to four.
#' @return MDRI Dataframe containing MDRI point estimates, CI lower and upper bounds and standard deviation of point estimates produced during bootstrapping. One row per functional form.
#' @return Plots A plot of Probability of testing recent over time for each functional form.
#' @return Models The fitted GLM models for each functional form.
#' @details The package contains long form documentation in the form of vignettes that cover the use of the main fucntions. Use browseVignettes(package="inctools") to access them.
#'
#' Expected data frame format: Before calling the function, please import your dataset into R environment.
#'
#' time_var: Time since infection; Note: this package does not assume any specific time unit. It is important to specify the
#' recency time cut-off 'T' and the time-based data exclusion rule (inclusion_time_threshold) in the same unit as the input times.
#' The estimated MDRI will be in this unit.
#'
#' Method: This function fits a function for probability of testing recent as a function of time to the supplied data using
#' binomial regression. This requires binary outcomes (recent/non-recent) coded as 1 for recent and 0 for non-recent test
#' resutls. Either a recency status variable must be specified, or a recency rule for determinging recency
#' status from a biomarker or set of biomarkers must be specified. Currently only independent biomarker
#' thresholds are supported (i.e. all biomarker criteria must be met in order for a specimen to be classified as recent).
#'
#' Functional forms currently supported for the binomial regression fitting procedure:
#' cloglog_linear, logit_cubic
#'
#' To be implemented in the near future: logit_spline
#'
#' logit_cubic: Fits a binomial regression to probability of testing recent with a logit link on a polynomial in t of the
#' third degree, where t is time since (detectable) infection.
#'
#' cloglog_linear: Fits a binomial regression to probability of testing recent with a log log link on log(t), where t is
#' time since (detectable) infection.
#'
#' recency_rule: binary_data - supply a binary variable with 1=recent and 0=non-recent in recency_vars.
#'
#' recency_rule:independent_thresholds: supply one threshold variable per biomarker in recency_vars and the relevant
#' thresholds, as well as whether a value below or above each threshold indicates recency in recency_params.
#'
#' recency_params expects a list of pairs of thresholds and thresholdtypes, with zero indicating a reading below the threshold
#'  implies recency and 1 that a reading above the threshold implies recency. (Note: two values, a threshold and a
#'  thresholdtype per variable must be specified in recency_params. For example, if you specify recency_vars =
#'  c('ODn','ViralLoad') you may specify recency_params = c(1.5,0,500,1), meaning that an ODn reading below 1.5 AND a
#'  viral load reasing above 500 indicates a recent result. Objects with missing values in its biomarker readings will be
#'  excluded from caculation.
#'
#' @examples
#' mdrical(data=excalibdata,
#'         subid_var = "SubjectID",
#'         time_var = "DaysSinceEDDI",
#'         recency_cutoff_time = 730.5,
#'         inclusion_time_threshold = 800,
#'         functional_forms = c("cloglog_linear"),
#'         recency_rule = "binary_data",
#'         recency_vars = "Recent",
#'         n_bootstraps = 100,
#'         alpha = 0.05,
#'         plot = TRUE)
#' @export
mdrical <- function(data = NULL, subid_var = NULL, time_var = NULL, functional_forms = c("cloglog_linear",
    "logit_cubic"), recency_cutoff_time = 730.5, inclusion_time_threshold = 800,
    recency_rule = "binary_data", recency_vars = NULL, recency_params = NULL,
    n_bootstraps = 100, alpha = 0.05, plot = TRUE, parallel = FALSE, cores = 4) {

  if (is.null(data)) {
    stop("No input data has been specified")
  }
  if (is.null(subid_var)) {
    stop("No subject identifier has been specified")
  }
  if (is.null(time_var)) {
    stop("No time variable has been specified")
  }

  if (is.null(recency_vars)) {
    stop("No recency variables have been specified")
  }

  if (!exists("data") || !is.data.frame(get("data"))) {
    stop("Specified data is not a dataframe or does not exist")
  }

  if (is.null(functional_forms)) {
    stop("Please select at least one functional form to apply to the data")
  }

  if (is.null(recency_rule)) {
    stop("Please specify a Recency Rule")
  }

  if (is.null(recency_vars)) {
    stop("Please specify at least one Recency Variable")
  }

  if (recency_rule == "binary_data") {
    if (length(recency_vars) > 1) {
      stop("Binary data should be specified in one recency (outcome) variable.")
    }
    if (!all(data$recency_vars == 0 | data$recency_vars == 1)) {
      stop("Input data is not binary")
    }
  }

  if (recency_rule == "independent_thresholds" & length(recency_vars) != 0.5 *
      length(recency_params)) {
    stop("The number of recency variables must match the number of recency paramaters.")
  }

  if (is.null(subid_var) | is.null(time_var)) {
    stop("Subject identifier and time variables must be specified.")
  }

  if (parallel == TRUE && Sys.info()["sysname"] == "Windows") {
    stop("Sorry, parallelisation of bootstrapping is not supported on Windows")
  }

  if (parallel == TRUE) {
    check_package("foreach")
    check_package("doMC")
  }

  # check that subject id, time and recency variables exist
  variables <- colnames(data)
  if (sum(variables == subid_var) != 1) {
    print(paste("There is no column", subid_var, "in the data frame."))
  }
  if (sum(variables == time_var) != 1) {
    print(paste("There is no column", time_var, "in the data frame."))
  }
  for (i in 1:length(recency_vars)) {
    if (sum(variables == recency_vars[i]) != 1) {
      print(paste("There is no column", recency_vars[i], "in the data frame."))
    }
  }

    ## Assign numeric subject ids, recency variables and recency status
    data <- process_data(data = data, subid_var = subid_var, time_var = time_var,
        recency_vars = recency_vars, inclusion_time_threshold = inclusion_time_threshold)
    data <- assign_recency_status(data = data, recency_params = recency_params, recency_rule = recency_rule)

    tolerance_glm2 = 1e-08
    maxit_glm2 = 50000
    tolerance_integral = 1e-08
    maxit_integral = 10000

    n_subjects <- max(data$sid)

    mdri_output <- data.frame(matrix(ncol = 4, nrow = 0))
    model_output <- list()
    if (plot == TRUE) {
        plot_output <- list()
    }

    for (i in 1:length(functional_forms)) {
        functional_form <- functional_forms[i]

        if (parallel == TRUE && n_bootstraps != 0) {
            boot_data <- data
            model <- fit_binomial_model(data = boot_data, functional_form = functional_form,
                tolerance = tolerance_glm2, maxit = maxit_glm2)
            parameters <- model$coefficients
            mdri <- integrate_for_mdri(parameters = parameters, recency_cutoff_time = recency_cutoff_time,
                functional_form = functional_form, tolerance = tolerance_integral,
                maxit = maxit_integral)
            mdris <- mdri
            model_output[[functional_forms[i]]] <- model
            if (plot == TRUE) {
                plot_parameters <- parameters
            }

            doMC::registerDoMC(cores)
            chosen_subjects <- vector(mode = "list", length = n_bootstraps)
            # set.seed(123)
            for (j in 1:n_bootstraps) {
                chosen_subjects[[j]] <- sample(1:n_subjects, n_subjects, replace = T)
            }
            mdris <- foreach::foreach(j = 1:n_bootstraps, .combine = rbind) %dopar%
                {
                  boot_data <- data[FALSE, ]
                  for (k in 1:n_subjects) {
                    boot_data <- rbind(boot_data, subset(data, sid == chosen_subjects[[j]][k]))
                  }
                  model <- fit_binomial_model(data = boot_data, functional_form = functional_form,
                    tolerance = tolerance_glm2, maxit = maxit_glm2)
                  parameters <- model$coefficients
                  mdri_iterate <- integrate_for_mdri(parameters = parameters, recency_cutoff_time = recency_cutoff_time,
                    functional_form = functional_form, tolerance = tolerance_integral,
                    maxit = maxit_integral)
                  return(mdri_iterate)
                }
        } else {
            # set.seed(123)
            for (j in 0:n_bootstraps) {
                chosen_subjects <- sample(1:n_subjects, n_subjects, replace = T)
                if (j != 0) {
                  boot_data <- data[FALSE, ]
                  for (k in 1:n_subjects) {
                    boot_data <- rbind(boot_data, subset(data, sid == chosen_subjects[k]))
                  }
                } else {
                  boot_data <- data
                }

                model <- fit_binomial_model(data = boot_data, functional_form = functional_form,
                  tolerance = tolerance_glm2, maxit = maxit_glm2)
                parameters <- model$coefficients
                mdri_iterate <- integrate_for_mdri(parameters = parameters, recency_cutoff_time = recency_cutoff_time,
                  functional_form = functional_form, tolerance = tolerance_integral,
                  maxit = maxit_integral)
                if (j == 0) {
                  mdri <- mdri_iterate
                  mdris <- mdri
                  model_output[[functional_forms[i]]] <- model
                  if (plot == TRUE) {
                    plot_parameters <- parameters
                  }
                } else {
                  mdris <- append(mdris, mdri_iterate)
                }
            }  # bootstraps
        }

        if (n_bootstraps == 0) {
            mdri_sd <- NA
            mdri_ci <- c(NA, NA)
            mdri_ff <- data.frame(round(mdri, 4), NA, NA, NA)
            mdri_output <- rbind(mdri_output, mdri_ff)
        } else {
            mdri_sd <- sd(mdris)
            mdri_ci <- quantile(mdris, probs = c(alpha/2, 1 - alpha/2))
            mdri_ff <- data.frame(round(mdri, 4), round(mdri_ci[1], 4), round(mdri_ci[2],
                4), round(mdri_sd, 4))
            mdri_output <- rbind(mdri_output, mdri_ff)
        }


        if (plot == TRUE) {
            plot_name <- functional_forms[i]
            plot_output[[plot_name]] <- plot_probability(functional_form = functional_form,
                parameters = plot_parameters, mdri = mdri, mdri_ci = mdri_ci, inclusion_time_threshold = inclusion_time_threshold,
                recency_cutoff_time = recency_cutoff_time)
        }

    }  # functional forms
    rownames(mdri_output) <- functional_forms
    colnames(mdri_output) <- c("PE", "CI_LB", "CI_UB", "SD")

    if (plot == TRUE) {
        output <- list(MDRI = mdri_output, Plots = plot_output, Models = model_output)
    } else {
        output <- list(MDRI = mdri_output, Models = model_output)
    }
    return(output)
}

# This is complicated - needs specification of the family of individual curves,
# function with parameters etc...  Estimate MDRI using a mixed effects binomial
# model...  mdri_ml_mixedbinomial <- function() { }



check_package <- function(package) {
    if (!require(package, character.only = TRUE)) {
        print(paste("Attempting to install dependency", package, sep = " "))
        install.packages(package, dependencies = TRUE)
        if (!require(package, character.only = TRUE)) {
            stop(paste("Package", package, "could not be automatically installed.",
                sep = " "))
        }
    }
}

process_data <- function(data = data, subid_var = subid_var, time_var = time_var,
    recency_vars = recency_vars, inclusion_time_threshold = inclusion_time_threshold) {
    names(data)[names(data) == subid_var] <- "sid"
    names(data)[names(data) == time_var] <- "time_since_eddi"
    temp_data <- data[, c("sid", "time_since_eddi")]
    for (i in 1:length(recency_vars)) {
        temp_data <- cbind(temp_data, data[, recency_vars[i]])
        colnames(temp_data)[2 + i] <- paste0("recency", i)
    }
    temp_data <- subset(temp_data, 0 < as.numeric(temp_data$time_since_eddi) & as.numeric(temp_data$time_since_eddi) <=
        inclusion_time_threshold)
    temp_data <- na.omit(temp_data)
    if (nrow(temp_data) < 1) {
        stop("Error: dataframe is empty after omitting rows with empty cells and applying time exclusion criterion")
    }
    data <- temp_data
    # replace non-numeric subject identifiers with unique numeric identifiers
    data$sid <- plyr::mapvalues(data$sid, from = unique(data$sid), to = seq(1:length(unique(data$sid))))
    # order by subject id and then time_since_eddi
    data$sid <- as.numeric(data$sid)
    data <- data[order(data$sid, data$time_since_eddi), ]
    return(data)
}

# Assign recency status to 0 and 1 using recency_vars and recency_params
assign_recency_status <- function(data = data, recency_params = recency_params, recency_rule = recency_rule) {

    switch(as.character(recency_rule), binary_data = {
        data$recency_status <- data$recency1
    }, independent_thresholds = {
        n_recvars <- length(recency_params)/2
        for (i in 1:n_recvars) {
            if (recency_params[2 * i] == 0) {
                data$recencytemp <- ifelse(data[, 2 + i] < recency_params[2 * i -
                  1], 1, 0)
            }
            if (recency_params[2 * i] == 1) {
                data$recencytemp <- ifelse(data[, 2 + i] > recency_params[2 * i -
                  1], 1, 0)
            }
            data <- plyr::rename(data, replace = c(recencytemp = paste0("recency_stat",
                i)))
        }
        data$recency_status <- ifelse(rowSums(data[(3 + n_recvars):ncol(data)]) >=
            n_recvars, 1, 0)
    })
    return(data)
}

fit_binomial_model <- function(data = data, functional_form = functional_form, tolerance, maxit) {
    data$time_since_eddi <- ifelse(data$time_since_eddi == 0, 1e-10, data$time_since_eddi)

    switch(as.character(functional_form), cloglog_linear = {
        fitted <- FALSE
        while (!fitted) {
            model <- glm2::glm2(formula = (1 - recency_status) ~ 1 + I(log(time_since_eddi)),
                family = binomial(link = "cloglog"), data = data, control = glm.control(epsilon = tolerance,
                  maxit = maxit, trace = FALSE))
            if (class(model)[1] == "try-error") {
                tolerance <- tolerance * 10
            } else {
                fitted <- TRUE
            }
        }
    }, logit_cubic = {
        fitted <- FALSE
        while (!fitted) {
            model <- glm2::glm2(formula = recency_status ~ 1 + I(time_since_eddi) +
                I(time_since_eddi^2) + I(time_since_eddi^3), family = binomial(link = "logit"),
                data = data, control = glm.control(epsilon = tolerance, maxit = maxit,
                  trace = FALSE))
            if (class(model)[1] == "try-error") {
                tolerance <- tolerance * 10
            } else {
                fitted <- TRUE
            }
        }
    })
    # coefficients <- model$coefficients
    return(model)
}

# The next two functions simply specify the model form for use in the integrator
functional_form_clogloglinear <- function(t, parameters) {
    exp(-exp(parameters[1] + (parameters[2]) * log(t)))
}

functional_form_logitcubic <- function(t, parameters) {
    1/(1 + exp(-(parameters[1] + parameters[2] * t + parameters[3] * t^2 + parameters[4] *
        t^3)))
}

# A function that integrates from 0 to T in order to obtain MDRI estimate
integrate_for_mdri <- function(parameters = parameters, recency_cutoff_time = recency_cutoff_time,
    functional_form = functional_form, tolerance, maxit) {

    if (is.nan(functional_form)) {
        stop("functional_form name required in order to evaluate functional form")
    }

    switch(as.character(functional_form), cloglog_linear = {
        answer <- try(cubature::adaptIntegrate(f = functional_form_clogloglinear,
            lowerLimit = 0, upperLimit = recency_cutoff_time, parameters = parameters,
            tol = tolerance, fDim = 1, maxEval = 0, absError = 0, doChecking = FALSE)$integral)
        if (class(answer) == "try-error") {
            cat("try-error", "\n")
            answer <- pracma::romberg(f = functional_form_clogloglinear, a = 0, b = recency_cutoff_time,
                parameters = parameters, tol = tolerance, maxit = maxit)$value
        }
    }, logit_cubic = {
        answer <- try(cubature::adaptIntegrate(f = functional_form_logitcubic, lowerLimit = 0,
            upperLimit = recency_cutoff_time, parameters = parameters, tol = tolerance,
            fDim = 1, maxEval = 0, absError = 0, doChecking = FALSE)$integral)
        if (class(answer) == "try-error") {
            cat("try-error", "\n")
            answer <- pracma::romberg(f = functional_form_logitcubic, a = 0, b = recency_cutoff_time,
                parameters = parameters, tol = tolerance, maxit = maxit)$value
        }
    })
    return(answer)
}

plot_probability <- function(functional_form = functional_form, parameters = parameters,
    mdri = mdri, inclusion_time_threshold = inclusion_time_threshold, recency_cutoff_time = recency_cutoff_time,
    mdri_ci = mdri_ci) {
    plot_time <- seq(from = 0, to = inclusion_time_threshold, by = 0.01)
    switch(as.character(functional_form), cloglog_linear = {
        plotdata <- data.frame(plot_time, exp(-exp(parameters[1] + (parameters[2]) *
            log(plot_time))))
        colnames(plotdata) <- c("time_since_eddi", "probability")
    }, logit_cubic = {
        plotdata <- data.frame(plot_time, 1/(1 + exp(-(parameters[1] + parameters[2] *
            plot_time + parameters[3] * plot_time^2 + parameters[4] * plot_time^3))))
        colnames(plotdata) <- c("time_since_eddi", "probability")
    })

    plotout <- ggplot2::ggplot() + ggplot2::geom_line(data = plotdata, ggplot2::aes(x = time_since_eddi,
        y = probability))
    plotout <- plotout + ggplot2::labs(x = "Time (since detectable infection)", y = "Probability of testing recent")
    plotout <- plotout + ggplot2::geom_vline(xintercept = mdri, colour = "blue")
    plotout <- plotout + ggplot2::geom_vline(xintercept = mdri_ci[1], colour = "blue",
        alpha = 0.7, linetype = "dashed")
    plotout <- plotout + ggplot2::geom_vline(xintercept = mdri_ci[2], colour = "blue",
        alpha = 0.7, linetype = "dashed")
    plotout <- plotout + ggplot2::annotate("text", label = "MDRI", x = mdri + 50,
        y = 0.95, colour = "blue")
    plotout <- plotout + ggplot2::geom_vline(xintercept = recency_cutoff_time, colour = "red")
    plotout <- plotout + ggplot2::annotate("text", label = "T", x = recency_cutoff_time +
        20, y = 0.95, colour = "red")
    plotout <- plotout + ggplot2::theme(panel.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = "dark grey"))
    plot_title <- paste0("Probability of testing recent over time (", functional_form,
        ")")
    plotout <- plotout + ggplot2::ggtitle(plot_title)

    return(plotout)
}
