library(inctools)
# No parallelisation 
system.time(print(mdrical(data=excalibdata,
        subid_var = "SubjectID",
        time_var = "DaysSinceEDDI",
        recency_cutoff_time = 730.5,
        inclusion_time_threshold = 800,
        functional_forms = c("logit_cubic"),
        recency_rule = "independent_thresholds",
        recency_vars = c("Result","VL"),
        recency_params = c(10,0,1000,1),
        n_bootstraps = 10000,
        alpha = 0.05,
        plot = FALSE,
        parallel = FALSE)))
# inctools 1.0.7: time: 982.879 SD:
# inctools 1.0.7.9100:

# With parallelisation
system.time(print(mdrical(data=excalibdata,
        subid_var = "SubjectID",
        time_var = "DaysSinceEDDI",
        recency_cutoff_time = 730.5,
        inclusion_time_threshold = 800,
        functional_forms = c("logit_cubic"),
        recency_rule = "independent_thresholds",
        recency_vars = c("Result","VL"),
        recency_params = c(10,0,1000,1),
        n_bootstraps = 10000,
        alpha = 0.05,
        plot = FALSE,
        parallel = TRUE,
        cores = 4)))
# inctools 1.0.7: time: 522.374 PE: 259.2234 SD: 14.5339
# inctools 1.0.7.9200: 48.630 PE: 259.2234 SD: 14.5723

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

system.time(mdri1 <- mdrical(data=excalibdata,
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
        plot = FALSE,
        parallel = TRUE,
        cores=4))
system.time(mdri2 <- mdrical(data=excalibdata,
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
                             plot = FALSE,
                             parallel = TRUE,
                             cores=4))

mdrical(data=excalibdata,
        subid_var = "SubjectID",
        time_var = "DaysSinceEDDI",
        recency_cutoff_time = 730.5,
        inclusion_time_threshold = 800,
        functional_forms = c("logit_cubic","cloglog_linear"),
        recency_rule = "independent_thresholds",
        recency_vars = c("Result","VL"),
        recency_params = c(10,0,1000,1),
        n_bootstraps = 100,
        alpha = 0.05,
        plot = FALSE,
        parallel = FALSE,
        cores=4,
        progress = FALSE)
