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

#' @importFrom magrittr "%>%"
#' @importFrom foreach "%dopar%"
#' @importFrom rlang .data
#' @importFrom rlang "!!"

# Function for resampling groups using dplyr
# inspired by drhagen
# https://github.com/tidyverse/dplyr/issues/361#issuecomment-243551042
# example use:
# iris %>% group_by(Species) %>% sample_frac_groups(1)
sample_frac_groups = function(tbl, size, replace = FALSE, weight=NULL) {
  # regroup when done
  grps <- tbl %>%
    dplyr::groups() %>%
    base::unlist() %>%
    base::as.character()
  # check length of groups non-zero
  keep <- tbl %>%
    dplyr::summarise() %>%
    dplyr::sample_frac(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping
  # variable
  tbl %>%
    dplyr::right_join(keep, by=grps) %>%
    dplyr::group_by(get(grps))
}

#' Estimate Mean Duration of Recent Infection (MDRI)
#'
#' Estimates MDRI (point estimate and confidence interval) using binomial
#' regression and a maximum likelihood approach
#'
#' @param data A data frame containing variables for subject identifier,
#' time (since detectable infection), and variables with biomarker readings or
#' recency status (to be specified in recency_vars)
#' @param functional_forms Select functional form/link function combinations
#' for fitting probability of testing recent as a function of time to data using
#' binomial regression
#' (see Details). Default=all supported functional forms.
#' @param subid_var The variable in the dataframe identifying subjects
#' @param time_var The variable in the dataframe indicating time between
#' 'time zero' (usually detectable infection) and biomarker measurement
# @param infwind_var The (optional) variable indicating the size of the
# seroconversion window, in the event time_var is relative to first positive
# test preceded by an infection window.  WE NEED A WINDOW PARAMATER AND WINDOW
# OPTIONS
# - e.g.  (1) uniform window width or (2) normal assumption, SDs
#' @param recency_rule Specified rule for defining recent/non-recent outcomes
#' from biomarker data (see Details)
#' @param recency_vars Variables to be used in determining recency outcomes
#' @param recency_params Vector of numeric parameters (e.g. thresholds) for
#' determining recency according to the relevant rule
#' @param recency_cutoff_time Recency time cut-off ('Big T'). Default = 730.5.
#' @param inclusion_time_threshold Data points beyond this time are excluded
#' from the calculation (in same unit as recency_cutoff_time, default = 800).
#' @param n_bootstraps Number of subject-level bootstrap resampling operations
#' for estimating confidence intervals, default = 10000.
#' @param random_seed Pass a random seed for reproducible bootstrapping.
#' Default is NULL.
#' @param alpha Confidence level, default=0.05.
# ADD OPTION TO GET FULL LIST OF MDRIs from the bootstrapping procedure or the
# shape of the distribution or something
#' @param plot Specifies whether a plot of the probability of testing recent
#' over time should be produced
#' @param parallel Set to TRUE in order to perform bootstrapping in parallel on
#' a multicore or multiprocessor system.
#' @param cores Set number of cores for parallel processing when parallel=TRUE.
#' This defaults to four.
#' @param output_bs_parms Return a matrix of the fitting parameters for each
#' bootstrap iteration.
#' @param debug Enable debugging mode (browser)
#' @return MDRI Dataframe containing MDRI point estimates, CI lower and upper
#' bounds and standard deviation of point estimates produced during
#' bootstrapping.
#' One row per functional form.
#' @return Plots A plot of Probability of testing recent over time for each
#' functional form.
#' @return Models The fitted GLM models for each functional form.
#' @details The package contains long form documentation in the form of
#' vignettes that cover the use of the main fucntions. Use
#' browseVignettes(package="inctools") to access them.
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
#'         n_bootstraps = 10,
#'         parallel = FALSE,
#'         alpha = 0.05,
#'         plot = TRUE)
#' @export
mdrical <- function(data = NULL,
                    subid_var = NULL,
                    time_var = NULL,
                    functional_forms = c("cloglog_linear", "logit_cubic"),
                    recency_cutoff_time = 730.5,
                    inclusion_time_threshold = 800,
                    recency_rule = "binary_data",
                    recency_vars = NULL,
                    recency_params = NULL,
                    n_bootstraps = 10000,
                    random_seed = NULL,
                    alpha = 0.05,
                    plot = TRUE,
                    parallel = ifelse(n_bootstraps == 0, FALSE, TRUE),
                    cores = parallel::detectCores(),
                    output_bs_parms = FALSE,
                    debug = FALSE) {

  if (debug) {browser()}

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

  if (!exists("data") || ( !is.data.frame(get("data")) & !tibble::is_tibble(get("data")) ) )  {
    stop("Specified data is not a dataframe or does not exist")
  }


  if (is.null(recency_rule)) {
    stop("Please specify a recency rule")
  }

  if (is.null(recency_vars)) {
    stop("Please specify at least one Recency Variable")
  }

    if (!(recency_rule %in% c("binary_data", "independent_thresholds"))) {
      stop("Please specify a valid recency rule")
    }

  if (recency_rule == "binary_data") {
    if (length(recency_vars) > 1) {
      stop("Binary data should have one recency variable")
    }
    # This line was broken. Should we have a similar check?
    # if (!all(data$recency_vars == 0 | data$recency_vars == 1)) {
    #   stop("Input data is not binary")
    # }
  }

  if (recency_rule == "independent_thresholds" & length(recency_vars) != 0.5 *
      length(recency_params)) {
    stop("The number of recency variables must match the number of recency paramaters")
  }

  if (is.null(functional_forms)) {
    stop("Please select at least one functional form to apply to the data")
  }
  
  for (ff in functional_forms) {
    if (!(ff %in% c("cloglog_linear", "logit_cubic"))) {
      stop("Please specify valid functional form(s)")
    }
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

  if (n_bootstraps < 0 | !is.wholenumber(n_bootstraps)) {
    stop("n_bootstraps must be a positive integer")
  }

  if (output_bs_parms & n_bootstraps == 0) {
    stop("Bootstrapped parameters can only be output if bootstrapping is performed")
  }

  if (parallel == TRUE & n_bootstraps == 0) {
    warning("Parallelisation only applicable if bootstrapping is performed")
    parallel <- FALSE
  }

  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  ## Assign numeric subject ids, recency variables and recency status
  data <- process_data(data = data,
                       subid_var = subid_var,
                       time_var = time_var,
                       recency_vars = recency_vars,
                       inclusion_time_threshold = inclusion_time_threshold,
                       debug = debug)
  data <- assign_recency_status(data = data,
                                recency_params = recency_params,
                                recency_rule = recency_rule,
                                debug = debug)

  tolerance_glm2 = 1e-08
  maxit_glm2 = 50000
  tolerance_integral = 1e-08
  maxit_integral = 10000

  #n_subjects <- max(data$sid)

  mdri_output <- tibble::tibble(
    FuncForm = character(),
    PE = numeric(),
    CI_LB = numeric(),
    CI_UB = numeric(),
    SE = numeric(),
    n_recent = numeric(),
    n_subjects = numeric(),
    n_observations = numeric(),
    .rows = 0
  ) 
  model_output <- list()

  if (output_bs_parms) {
    bs_parms_output <- list()
  }

  if (plot == TRUE) {
    plot_output <- list()
  }

  for (i in 1:length(functional_forms)) {
    functional_form <- functional_forms[i]
    #print(paste("Computing MDRI using functional form",functional_form))

    if (parallel == TRUE && n_bootstraps > 0) {
      model <- fit_binomial_model(data = data,
                                  functional_form = functional_form,
                                  tolerance = tolerance_glm2,
                                  maxit = maxit_glm2)
      parameters <- model$coefficients
      mdri <- integrate_for_mdri(parameters = parameters,
                                 recency_cutoff_time = recency_cutoff_time,
                                 functional_form = functional_form,
                                 tolerance = tolerance_integral,
                                 maxit = maxit_integral)
      mdris <- mdri
      model_output[[functional_forms[i]]] <- model
      if (plot == TRUE) {
        plot_parameters <- parameters
      }
      ## WORKAROUND: https://github.com/rstudio/rstudio/issues/6692
      ## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4
      if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
          Sys.info()["sysname"] == "Darwin" && R.Version()$major == "4") {
        cluster <- parallel::makeCluster(cores, outfile="", setup_strategy = "sequential")
      } else {
        cluster <- parallel::makeCluster(cores, outfile="")
      }
      
      doParallel::registerDoParallel(cluster)
      if (foreach::getDoParWorkers() != cores) {
        stop("Failed to initialise parallel worker threads.")
      }
      pb <- utils::txtProgressBar(min = 1, max = n_bootstraps, style = 3)

      # Group data for bootstrapping purposes
      data_grouped <- data %>%
        dplyr::group_by(.data$sid)

      if (!output_bs_parms) {
        mdris <- foreach::foreach(j = 1:n_bootstraps, .combine = c,
                                  #.options.snow = opts,
                                  .inorder = FALSE #,
                                  #.packages = "inctools"
        ) %dopar%
        {
          if (!is.null(random_seed)) {
            set.seed(j * random_seed)
          }
          boot_data <- data_grouped %>%
            sample_frac_groups(1, replace = TRUE) %>%
            dplyr::ungroup()
          model <- fit_binomial_model(data = boot_data,
                                      functional_form = functional_form,
                                      tolerance = tolerance_glm2,
                                      maxit = maxit_glm2)
          parameters <- model$coefficients
          mdri_iterate <- integrate_for_mdri(parameters = parameters,
                                             recency_cutoff_time = recency_cutoff_time,
                                             functional_form = functional_form,
                                             tolerance = tolerance_integral,
                                             maxit = maxit_integral)
          if (n_bootstraps > 0) {utils::setTxtProgressBar(pb, j)}
          return(mdri_iterate)
        }
        close(pb)
        parallel::stopCluster(cluster)

      } else if(output_bs_parms) {
        mdris_and_params <- foreach::foreach(j = 1:n_bootstraps, .combine = dplyr::bind_rows,
                                             #.options.snow = opts,
                                             .inorder = FALSE #,
                                             #.packages = "inctools"
        ) %dopar%
        {
          if (!is.null(random_seed)) {
            set.seed(j * random_seed)
          }
          boot_data <- data_grouped %>%
            sample_frac_groups(1, replace = TRUE) %>%
            dplyr::ungroup()
          model <- fit_binomial_model(data = boot_data,
                                      functional_form = functional_form,
                                      tolerance = tolerance_glm2,
                                      maxit = maxit_glm2)
          parameters <- model$coefficients

          if(length(parameters) == 4) {
            names(parameters) <- c("beta0","beta1","beta2","beta3")
          } else if (length(parameters) == 2) {
            names(parameters) <- c("beta0","beta1")
          }

          mdri_iterate <- integrate_for_mdri(parameters = parameters,
                                             recency_cutoff_time = recency_cutoff_time,
                                             functional_form = functional_form,
                                             tolerance = tolerance_integral,
                                             maxit = maxit_integral)
          if (n_bootstraps > 0) {utils::setTxtProgressBar(pb, j)}

          mdri_and_params_iterate <- tibble::tibble("MDRI" = mdri_iterate) %>%
            dplyr::bind_cols(tibble::as_tibble(t(parameters)))

          return(mdri_and_params_iterate)
        }
        close(pb)
        parallel::stopCluster(cluster)

        mdris <- as.vector(mdris_and_params$MDRI)

        bs_params <- mdris_and_params %>%
          dplyr::select(-"MDRI")

        bs_parms_output[[functional_forms[i]]] <- bs_params
      }

    } else if(!parallel) {
      if (n_bootstraps > 0) {pb <- utils::txtProgressBar(min = 1, max = n_bootstraps, style = 3)}

      # Group data for bootstrapping purposes
      data_grouped <- data %>%
        dplyr::group_by(.data$sid)

      for (j in 0:n_bootstraps) {
        if (j != 0) {
          boot_data <- data_grouped %>%
            sample_frac_groups(1, replace = TRUE) %>%
            dplyr::ungroup()
        } else {
          boot_data <- data
        }

        model <- fit_binomial_model(data = boot_data,
                                    functional_form = functional_form,
                                    tolerance = tolerance_glm2,
                                    maxit = maxit_glm2)
        parameters <- model$coefficients
        if(length(parameters) == 4) {
          names(parameters) <- c("beta0","beta1","beta2","beta3")
        } else if (length(parameters) == 2) {
          names(parameters) <- c("beta0","beta1")
        }

        mdri_iterate <- integrate_for_mdri(parameters = parameters,
                                           recency_cutoff_time = recency_cutoff_time,
                                           functional_form = functional_form,
                                           tolerance = tolerance_integral,
                                           maxit = maxit_integral)
        if (j == 0) {
          mdri <- mdri_iterate
          mdris <- mdri
          model_output[[functional_forms[i]]] <- model
          if (plot == TRUE) {
            plot_parameters <- parameters
          }
          if(output_bs_parms) {
            bs_params <- tibble::as_tibble(t(parameters))[NULL,]
            }
        } else if(j > 0) {
          mdris <- append(mdris, mdri_iterate)
          if(output_bs_parms) {
          bs_params <- bs_params %>%
            dplyr::bind_rows(tibble::as_tibble(t(parameters)))
          bs_parms_output[[functional_forms[i]]] <- bs_params
          }
          utils::setTxtProgressBar(pb, j)
        }
      }  # bootstraps
      if (n_bootstraps > 0) {close(pb)}
    }

    if (n_bootstraps == 0) {
      mdri_sd <- NA
      mdri_ci <- c(NA, NA)
      mdri_output <- mdri_output %>%
        dplyr::bind_rows(
          tibble::tibble(
            FuncForm = functional_form,
            PE = mdri,
            CI_LB = mdri_ci[1],
            CI_UB = mdri_ci[2],
            SE = mdri_sd,
            n_recent = sum(data$recency_status),
            n_subjects = length(unique(data$sid)),
            n_observations = nrow(data)
          ) 
        )
    } else {
      mdri_sd <- stats::sd(mdris)
      mdri_ci <- stats::quantile(mdris, probs = c(alpha/2, 1 - alpha/2))
      mdri_output <- mdri_output %>%
        dplyr::bind_rows(
          tibble::tibble(
            FuncForm = functional_form,
            PE = mdri,
            CI_LB = mdri_ci[1],
            CI_UB = mdri_ci[2],
            SE = mdri_sd,
            n_recent = sum(data$recency_status),
            n_subjects = length(unique(data$sid)),
            n_observations = nrow(data)
          ) 
        )
    }


    if (plot == TRUE) {
      plot_name <- functional_forms[i]
      plot_output[[plot_name]] <- plot_probability(functional_form = functional_form,
                                                   parameters = plot_parameters, mdri = mdri, mdri_ci = mdri_ci, inclusion_time_threshold = inclusion_time_threshold,
                                                   recency_cutoff_time = recency_cutoff_time)
    }

  }  # functional forms

  if (!plot) {plot_output <- NULL}
  if (!output_bs_parms) {bs_parms_output <- NULL}

  output <- list(MDRI = mdri_output, Models = model_output, Plots = plot_output, BSparms = bs_parms_output)

  # Does this affect the calling environment?
  options(pillar.sigfig = 6)
  return(output)
}

process_data <- function(data = data,
                         subid_var = subid_var,
                         time_var = time_var,
                         recency_vars = recency_vars,
                         inclusion_time_threshold = inclusion_time_threshold,
                         debug = FALSE) {

  if (debug) {browser()}
  
  recency_vars_newnames <- paste0("recency", 1:length(recency_vars))
  data <- data %>%
    dplyr::rename(sid = !!subid_var,
           time_since_eddi = !!time_var) %>%
    dplyr::select(.data$sid, .data$time_since_eddi, recency_vars) %>%
    dplyr::rename_at(recency_vars, function(x) recency_vars_newnames) %>%
    dplyr::filter(.data$time_since_eddi > 0, .data$time_since_eddi <= inclusion_time_threshold) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(time_since_eddi = as.numeric(as.character(.data$time_since_eddi)),
                  sid = plyr::mapvalues(.data$sid, unique(.data$sid), seq(1:length(unique(.data$sid))))) %>% # Is there a better way than using plyr?
    dplyr::arrange(.data$sid, .data$time_since_eddi)
  
  if (nrow(data) < 1) {
    stop("Error: dataframe is empty after omitting rows with empty cells and applying time exclusion criterion")
  }
  
  return(data)
}

# Assign recency status to 0 and 1 using recency_vars and recency_params
assign_recency_status <- function(data = data,
                                  recency_params = recency_params,
                                  recency_rule = recency_rule,
                                  debug = FALSE) {
  if(debug) {browser()}

  
  if (recency_rule == "binary_data") {
    
    data <- dplyr::rename(data, recency_status = .data$recency1)
  
  } else if (recency_rule == "independent_thresholds") {
   
    n_recentvars <- length(recency_params)/2
    recencyvars <- paste0("recency", 1:n_recentvars)
    statusvars <- paste0("recency_stat", 1:n_recentvars)
    data <- data %>%
      dplyr::rename_at(recencyvars, function(x) statusvars)
   
    for (i in 1:n_recentvars) {
      if (recency_params[2 * i] == 0) {
        
        data[,2+i] <- dplyr::case_when(
            data[,2+i] < recency_params[2 * i - 1] ~ 1,
            data[,2+i] >= recency_params[2 * i - 1] ~ 0
            )
        
      }
      
      if (recency_params[2 * i] == 1) {
        
        data[,2+i] <- dplyr::case_when(
          data[,2+i] > recency_params[2 * i - 1] ~ 1,
          data[,2+i] <= recency_params[2 * i - 1] ~ 0
        )
        
      }
    }
    
    data <- data %>%
      dplyr::mutate(recency_sum = rowSums(data[,3:(2+n_recentvars)]),
             recency_status = dplyr::case_when(
               .data$recency_sum < n_recentvars ~ 0,
               .data$recency_sum >= n_recentvars ~ 1,
             ))
  } else {
    
    stop("Error: Recency rule is not `binary_data` or `independent_thresholds`.")
    
  }
  
  return(data)
}

fit_binomial_model <- function(data = data,
                               functional_form = functional_form,
                               tolerance = tolerance,
                               maxit = maxit) {
  data$time_since_eddi <- ifelse(data$time_since_eddi == 0,
                                 1e-10,
                                 data$time_since_eddi)

  switch(as.character(functional_form), cloglog_linear = {
    fitted <- FALSE
    while (!fitted) {
      suppressWarnings(model <- glm2::glm2(formula = recency_status ~ 1 +
                                             I(log(time_since_eddi)),
                                           family = stats::binomial(link = "cloglog"),
                                           data = data,
                                           control = stats::glm.control(epsilon = tolerance,
                                                                        maxit = maxit,
                                                                        trace = FALSE)))
      if (class(model)[1] == "try-error") {
        tolerance <- tolerance * 10
      } else {
        fitted <- TRUE
      }
    }
  }, logit_cubic = {
    fitted <- FALSE
    while (!fitted) {
      suppressWarnings(model <- glm2::glm2(formula = recency_status ~ 1 + I(time_since_eddi) +
                                             I(time_since_eddi^2) + I(time_since_eddi^3),
                                           family = stats::binomial(link = "logit"),
                                           data = data,
                                           control = stats::glm.control(epsilon = tolerance,
                                                                        maxit = maxit,
                                                                        trace = FALSE)))
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
  1 - exp(-exp(parameters[1] + (parameters[2]) * log(t)))
}

functional_form_logitcubic <- function(t, parameters) {
  1/(1 + exp(-(parameters[1] + parameters[2] * t + parameters[3] * t^2 +
                 parameters[4] * t^3)))
}

integrate_for_mdri <- function(parameters = parameters,
                               recency_cutoff_time = recency_cutoff_time,
                               functional_form = functional_form,
                               tolerance,
                               maxit) {

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

plot_probability <- function(functional_form = functional_form,
                             parameters = parameters,
                             mdri = mdri,
                             inclusion_time_threshold = inclusion_time_threshold,
                             recency_cutoff_time = recency_cutoff_time,
                             mdri_ci = mdri_ci) {
  plot_time <- seq(from = 0, to = inclusion_time_threshold, by = 0.01)
  switch(as.character(functional_form), cloglog_linear = {
    plotdata <- tibble::tibble(plot_time, functional_form_clogloglinear(t = plot_time, parameters = parameters))
    colnames(plotdata) <- c("time_since_eddi", "probability")
  }, logit_cubic = {
    plotdata <- tibble::tibble(plot_time, functional_form_logitcubic(t = plot_time, parameters = parameters))
    colnames(plotdata) <- c("time_since_eddi", "probability")
  })

  plotout <- ggplot2::ggplot() +
    ggplot2::geom_line(data = plotdata,
                       ggplot2::aes(x = .data$time_since_eddi,
                                    y = .data$probability)) +
    ggplot2::labs(x = "Time (since detectable infection)",
                  y = "Probability of testing recent") +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    ggplot2::geom_vline(xintercept = mdri, colour = "blue")

  if (!is.na(mdri_ci[1]) & !is.na(mdri_ci[2]) & !is.null(mdri_ci[1]) &
      !is.null(mdri_ci[2])) {
    plotout <- plotout + ggplot2::geom_vline(xintercept = mdri_ci[1],
                                             colour = "blue", alpha = 0.7,
                                             linetype = "dashed") +
      ggplot2::geom_vline(xintercept = mdri_ci[2],
                          colour = "blue", alpha = 0.7,
                          linetype = "dashed")
  }
  plotout <- plotout + ggplot2::annotate("text", label = "MDRI", x = mdri + 50,
                                         y = 0.95, colour = "blue") +
    ggplot2::geom_vline(xintercept = recency_cutoff_time, colour = "red") +
    ggplot2::annotate("text", label = "T", x = recency_cutoff_time +
                        20, y = 0.95, colour = "red") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "dark grey")) +
    ggplot2::ggtitle(paste0("Probability of testing recent over time (", functional_form,
                            ")"))

  return(plotout)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
