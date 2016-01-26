# Install devtools
install.packages("devtools")

# Install dependencies manually (shouldn't be necessary ordinarily)
install.packages(c("plyr","glm2","cubature","MASS","ggplot2","pracma"))
# If you are not on Windows and wish to perform bootstrapping in parallel, uncomment the following line
# install.packages(c("foreach","doMC"))

# Install package using devtools
# Note: Change the working directory in the following line to the local location of "ritcalib"
setwd("~/Development/ABIE/abie")
devtools::install(dependencies=TRUE)

# Load the package in order to use it
library(ritcalib)

setwd("~/Development/Recent-Infection-Test-Calibration/ritcalib")
exampledata <- read.csv("data/exampledata.csv")

# In the following example, subject (panel) identifier and time since estimated date of infection
# is provided. Then a list of "functional forms" - this parameter can be left out since the default is
# to use both.
# The next two paramters are the recency cutoff time aka "Big T" (we usually suggest 2 years) and the
# data exclusion rule, i.e. all data points beyond this time are ignored in the fitting procedure.
# The "recency rule" is to tell the function how to distinguish recent from non-recent results. You
# can classify the data yourself and use "binary_data" as the rule, and then specify the variable in
# your data frame
# We are also specifying a vector of variables and a vector of paramaters to define recency
# In this case we are using the assay result and the viral load. The paramaters 10,0,1000,1 mean
# that recency is defined as a biomarker reading below 10 and a viral load reading above 1000
# alpha and n_bootstraps relate to obtaining confidence intervals using subject-level bootstrapping.

mdri_ml_binomial(data=exampledata,
                 subid_var = "SubjectID",
                 time_var = "DaysSinceEDDI",
                 recency_cutoff_time = 730.5,
                 inclusion_time_threshold = 800,
                 recency_rule = "independent_thresholds",
                 recency_vars = c("Result","VL"),
                 recency_params = c(10,0,1000,1),
                 n_bootstraps = 100,
                 alpha = 0.05,
                 plot = TRUE)

# As above, but parellelise the bootstrapping. In this case, split the job over four cores.

mdri_ml_binomial(data=exampledata,
                 subid_var = "SubjectID",
                 time_var = "DaysSinceEDDI",
                 recency_cutoff_time = 730.5,
                 inclusion_time_threshold = 800,
                 recency_rule = "independent_thresholds",
                 recency_vars = c("Result","VL"),
                 recency_params = c(10,0,1000,1),
                 n_bootstraps = 100,
                 alpha = 0.05,
                 plot = TRUE,
                 parallel = TRUE,
                 cores=4)

# This example calculates a false-recent rate, treating the data at subject level
frr_binomial(data=exampledata,
             subid_var = "SubjectID",
             time_var = "DaysSinceEDDI",
             recency_cutoff_time = 730.5,
             recency_rule = "independent_thresholds",
             recency_vars = c("Result","VL"),
             recency_params = c(10,0,1000,1),
             alpha = 0.05)
