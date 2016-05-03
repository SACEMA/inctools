## ------------------------------------------------------------------------
library(abie)

## ------------------------------------------------------------------------
exampledata <- read.csv("../data/exampledata_testcalibration.csv")

## ---- fig.width=6.5, fig.height=5, fig.align="center"--------------------
mdrical(data=exampledata,
                 subid_var = "SubjectID",
                 time_var = "DaysSinceEDDI",
                 recency_cutoff_time = 730.5,
                 inclusion_time_threshold = 800,
                 functional_forms = c("cloglog_linear"),
                 recency_rule = "independent_thresholds",
                 recency_vars = c("Result","VL"),
                 recency_params = c(10,0,1000,1),
                 n_bootstraps = 100,
                 alpha = 0.05,
                 plot = TRUE)

## ---- fig.width=6.5, fig.height=5, fig.align="center"--------------------
mdrical(data=exampledata,
                 subid_var = "SubjectID",
                 time_var = "DaysSinceEDDI",
                 recency_cutoff_time = 730.5,
                 inclusion_time_threshold = 800,
                 functional_forms = c("logit_cubic"),
                 recency_rule = "independent_thresholds",
                 recency_vars = c("Result","VL"),
                 recency_params = c(10,0,1000,1),
                 n_bootstraps = 100,
                 alpha = 0.05,
                 plot = TRUE)

## ---- fig.width=7, fig.height=5, fig.align="center"----------------------
mdrical(data=exampledata,
                 subid_var = "SubjectID",
                 time_var = "DaysSinceEDDI",
                 recency_cutoff_time = 730.5,
                 inclusion_time_threshold = 800,
                 functional_forms = c("logit_cubic","cloglog_linear"),
                 recency_rule = "independent_thresholds",
                 recency_vars = c("Result","VL"),
                 recency_params = c(10,0,1000,1),
                 n_bootstraps = 1000,
                 alpha = 0.05,
                 plot = TRUE,
                 parallel = TRUE,
                 cores=4)

## ------------------------------------------------------------------------
frrcal(data=exampledata,
             subid_var = "SubjectID",
             time_var = "DaysSinceEDDI",
             recency_cutoff_time = 730.5,
             recency_rule = "independent_thresholds",
             recency_vars = c("Result","VL"),
             recency_params = c(10,0,1000,1),
             alpha = 0.05)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

