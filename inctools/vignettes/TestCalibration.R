## ------------------------------------------------------------------------
library(inctools)

## ---- fig.width=6.5, fig.height=5, fig.align="center", fig.show='hold'----
mdrical(data=excalibdata,
                 subid_var = "SubjectID",
                 time_var = "DaysSinceEDDI",
                 recency_cutoff_time = 730.5,
                 inclusion_time_threshold = 800,
                 functional_forms = c("cloglog_linear"),
                 recency_rule = "binary_data",
                 recency_vars = "Recent",
                 n_bootstraps = 1000,
                 alpha = 0.05,
                 plot = TRUE)

## ---- fig.width=6.5, fig.height=5, fig.align="center", fig.show='hold'----
mdrical(data=excalibdata,
                 subid_var = "SubjectID",
                 time_var = "DaysSinceEDDI",
                 recency_cutoff_time = 730.5,
                 inclusion_time_threshold = 800,
                 functional_forms = c("logit_cubic"),
                 recency_rule = "independent_thresholds",
                 recency_vars = c("Result","VL"),
                 recency_params = c(10,0,1000,1),
                 n_bootstraps = 1000,
                 alpha = 0.05,
                 plot = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  mdrical(data=excalibdata,
#                   subid_var = "SubjectID",
#                   time_var = "DaysSinceEDDI",
#                   recency_cutoff_time = 730.5,
#                   inclusion_time_threshold = 800,
#                   functional_forms = c("logit_cubic","cloglog_linear"),
#                   recency_rule = "independent_thresholds",
#                   recency_vars = c("Result","VL"),
#                   recency_params = c(10,0,1000,1),
#                   n_bootstraps = 10000,
#                   alpha = 0.05,
#                   plot = TRUE,
#                   parallel = TRUE,
#                   cores=4)

## ------------------------------------------------------------------------
frrcal(data=excalibdata,
             subid_var = "SubjectID",
             time_var = "DaysSinceEDDI",
             recency_cutoff_time = 730.5,
             recency_rule = "independent_thresholds",
             recency_vars = c("Result","VL"),
             recency_params = c(10,0,1000,1),
             alpha = 0.05)

