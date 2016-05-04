## ------------------------------------------------------------------------
library(inctools)

## ------------------------------------------------------------------------
exampledata <- read.csv("../data/exampledata_testcalibration.csv")

## ---- fig.width=6.5, fig.height=5, fig.align="center", fig.show='hold'----
mdrical(data=exampledata,
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

