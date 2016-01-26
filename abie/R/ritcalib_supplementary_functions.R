# Assay-Based Incidence Estimation
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
#
check_package <- function(package) {
  if (!require(package,character.only = TRUE)) {
    print(paste("Attempting to install dependency",package,sep=" "))
    install.packages(package,dep=TRUE)
    if (!require(package,character.only = TRUE)) {
      stop(paste("Package",package,"could not be automatically installed.",sep=" "))
    }
  }
}

process_data<-function(data=data, subid_var=subid_var, time_var=time_var, recency_vars=recency_vars, inclusion_time_threshold=inclusion_time_threshold){
  names(data)[names(data)==subid_var] <- "sid"
  names(data)[names(data)==time_var]   <- "time_since_eddi"
  temp_data <- data[ ,c("sid","time_since_eddi")]
  for (i in 1:length(recency_vars)) {
    temp_data <- cbind(temp_data, data[ , recency_vars[i]])
    colnames(temp_data)[2+i] <- paste0("recency", i)
  }
  temp_data <- subset(temp_data,0<as.numeric(temp_data$time_since_eddi) &
                        as.numeric(temp_data$time_since_eddi)<=inclusion_time_threshold)
  temp_data <- na.omit(temp_data)
  if (nrow(temp_data)<1)
    {stop("Error: dataframe is empty after omitting rows with empty cells and applying time exclusion criterion")}
  data <- temp_data
  # replace non-numeric subject identifiers with unique numeric identifiers
  data$sid <- plyr::mapvalues(data$sid, from = unique(data$sid), to = seq(1:length(unique(data$sid))))
  # order by subject id and then time_since_eddi
  data$sid <- as.numeric(data$sid)
  data <- data[order(data$sid,data$time_since_eddi), ]
  return(data)
}

#Assign recency status to 0 and 1 using recency_vars and recency_params
assign_recency_status <- function(data=data, recency_params=recency_params, recency_rule=recency_rule) {

  switch(as.character(recency_rule),
         binary_data = {
           data$recency_status <- data$recency1
         },
         independent_thresholds = {
           n_recvars <- length(recency_params)/2
           for (i in 1:n_recvars) {
             if (recency_params[2*i] == 0) {data$recencytemp <- ifelse (data[, 2+i] < recency_params[2*i-1],1,0)}
             if (recency_params[2*i] == 1) {data$recencytemp <- ifelse (data[, 2+i] > recency_params[2*i-1],1,0)}
             data <- plyr::rename(data,replace=c("recencytemp" = paste0("recency_stat",i)))
           }
           data$recency_status <- ifelse(rowSums(data[(3+n_recvars):ncol(data)])>=n_recvars,1,0)
         })
  return(data)
}

fit_binomial_model <- function(data=data, functional_form=functional_form,tolerance=tolerance_glm2, maxit=maxit_glm2) {
  data$time_since_eddi <- ifelse(data$time_since_eddi == 0, 1e-10, data$time_since_eddi)

  switch(as.character(functional_form),
         loglog_linear = {
           fitted <- FALSE
           while (!fitted) {
             model <- glm2::glm2(formula = (1-recency_status) ~ 1 + I(log(time_since_eddi)),
                           family = binomial(link = "cloglog"), data = data,
                           control = glm.control(epsilon = tolerance, maxit = maxit, trace = FALSE))
             if(class(model)[1] == "try-error"){
               tolerance <- tolerance*10
             }else{
               fitted <- TRUE
             }}
         },
         logit_cubic = {
           fitted <- FALSE
           while(!fitted){
             model <- glm2::glm2(formula = recency_status ~ 1 + I(time_since_eddi) + I(time_since_eddi^2) + I(time_since_eddi^3),
                           family = binomial(link = "logit"), data = data,
                           control = glm.control(epsilon = tolerance, maxit = maxit, trace = FALSE))
             if(class(model)[1] == "try-error"){
               tolerance <- tolerance*10
             }else{
               fitted <- TRUE
             }}
         })
  #coefficients <- model$coefficients
  return(model)
}

# The next two functions simply specify the model form for use in the integrator
functional_form_clogloglinear <- function(t, parameters) {
    exp(-exp(parameters[1] + (parameters[2])*log(t)))
}

functional_form_logitcubic <- function(t, parameters) {
    1/(1+exp(-(parameters[1]+parameters[2]*t+parameters[3]*t^2+parameters[4]*t^3)))
}

# A function that integrates from 0 to T in order to obtain MDRI estimate
integrate_for_mdri<-function(parameters=parameters, recency_cutoff_time = recency_cutoff_time, functional_form=functional_form, tolerance=tolerance_integral, maxit=maxit_integral){

  if (is.nan(functional_form)) {stop("functional_form name required in order to evaluate functional form")}

  switch(as.character(functional_form),
         loglog_linear = {
           answer <- try(cubature::adaptIntegrate(f = functional_form_clogloglinear, lowerLimit = 0,
                                                  upperLimit = recency_cutoff_time, parameters=parameters,
                                                  tol = tolerance, fDim = 1, maxEval = 0, absError = 0,
                                                  doChecking = FALSE)$integral)
           if(class(answer)=="try-error"){
             cat("try-error","\n")
             answer <- pracma::romberg(f = functional_form_clogloglinear, a = 0, b = recency_cutoff_time,
                                       parameters=parameters, tol = tolerance, maxit = maxit)$value}
         },
         logit_cubic = {
           answer <- try(cubature::adaptIntegrate(f = functional_form_logitcubic, lowerLimit = 0,
                                                  upperLimit = recency_cutoff_time, parameters=parameters,
                                                  tol = tolerance, fDim = 1, maxEval = 0, absError = 0,
                                                  doChecking = FALSE)$integral)
           if(class(answer)=="try-error"){
             cat("try-error","\n")
             answer <- pracma::romberg(f = functional_form_logitcubic, a = 0, b = recency_cutoff_time,
                                       parameters=parameters, tol = tolerance, maxit = maxit)$value}
         })
  return(answer)
}

plot_probability<-function(functional_form=functional_form, parameters=parameters, mdri=mdri,
                           inclusion_time_threshold=inclusion_time_threshold, recency_cutoff_time=recency_cutoff_time,
                           mdri_ci=mdri_ci) {
  plot_time <- seq(from=0,to=inclusion_time_threshold,by=0.01)
  switch(as.character(functional_form),
         loglog_linear = {
           plotdata<-data.frame(plot_time,
                                exp(-exp(parameters[1] + (parameters[2])*log(plot_time))))
           colnames(plotdata)<-c("time_since_eddi", "probability")
         },
         logit_cubic = {
           plotdata<-data.frame(plot_time,
                                1/(1+exp(-(parameters[1]+parameters[2]*plot_time+parameters[3]*plot_time^2+parameters[4]*plot_time^3))))
           colnames(plotdata)<-c("time_since_eddi", "probability")
         })

  plotout <- ggplot2::ggplot() + ggplot2::geom_line(data=plotdata, ggplot2::aes(x=time_since_eddi, y=probability))
  plotout <- plotout + ggplot2::labs(x="Time (since detectable infection)",y="Probability of testing recent")
  plotout <- plotout + ggplot2::geom_vline(xintercept=mdri, colour="blue")
  plotout <- plotout + ggplot2::geom_vline(xintercept=mdri_ci[1], colour="blue", alpha=0.7, linetype="dashed")
  plotout <- plotout + ggplot2::geom_vline(xintercept=mdri_ci[2], colour="blue", alpha=0.7, linetype="dashed")
  plotout <- plotout + ggplot2::annotate("text",label="MDRI", x=mdri+50 ,y=0.95, colour="blue")
  plotout <- plotout + ggplot2::geom_vline(xintercept=recency_cutoff_time, colour="red")
  plotout <- plotout + ggplot2::annotate("text",label="T", x=recency_cutoff_time+20 , y=0.95, colour="red")
  plotout <- plotout + ggplot2::theme(panel.background=ggplot2::element_blank(),panel.grid.major=ggplot2::element_line(colour="dark grey"))
  plot_title <- paste0("Probability of testing recent over time (",functional_form,")")
  plotout <- plotout + ggplot2::ggtitle(plot_title)

  return(plotout)}
