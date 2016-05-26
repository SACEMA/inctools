
# Try the following with versions of the package using doMC and doSNOW

system.time(mdrical(data=excalibdata,
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
                    plot = FALSE,
                    parallel = FALSE))

system.time(mdrical(data=excalibdata,
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
        plot = FALSE,
        parallel = TRUE,
        cores=4))


