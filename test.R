n_cores <- 4
setwd("~/dev/incidence-kzn/")

library(dplyr); library(tidyr);

lag <- read.csv("data/20170116-EP-LAgSedia.csv", stringsAsFactors = FALSE)
lag <- filter(lag, result_field=="final_result" & visit_hivstatus=="P")
lag <- rename(lag, LAg_ODn=result_value)
lag <- separate(lag, specimen_label, c("speclabel_root","speclabel_suffix"), sep = "-", remove = FALSE)
biorad <- read.csv("data/20160927_EP_BioradAvidityCDC.csv", stringsAsFactors = FALSE)
biorad <- filter(biorad, result_field=="final_result" & visit_hivstatus=="P")
biorad <- rename(biorad, BioRadAI=result_value)
biorad <- separate(biorad, specimen_label, c("speclabel_root","speclabel_suffix"), sep = "-", remove = FALSE)
biorad <- select(biorad, speclabel_root, BioRadAI)
lagbr <- left_join(lag, biorad)

lagbr <- lagbr %>%
  mutate(days_since_lp_ddi = ifelse(is.na(days_since_lp_ddi) &
                                      !is.na(days_since_cohort_entry) &
                                      cohort_entry_hiv_status == "P", days_since_cohort_entry, days_since_lp_ddi))

untreated_c <- dplyr::filter(lagbr, treatment_naive=="True" & scopevisit_ec=="False" & subtype == "C")
untreated_all <- dplyr::filter(lagbr, treatment_naive=="True" & scopevisit_ec=="False")
treated_all <- dplyr::filter(lagbr, on_treatment=="True" & first_treatment=="True" & days_since_current_art>=90 & scopevisit_ec=="False")
treated_c <- dplyr::filter(lagbr, on_treatment=="True" & first_treatment=="True" & days_since_current_art>=90 & scopevisit_ec=="False" & subtype == "C")



mdri <- mdrical(data = untreated_c,
                subid_var = "subject_id",
                time_var = "days_since_eddi",
                functional_forms = "logit_cubic",
                recency_cutoff_time = 730.5,
                recency_rule = "independent_thresholds",
                recency_vars = c("LAg_ODn","viral_load"),
                recency_params = c(1.5,0,75,1),
                n_bootstraps = 10000,
                plot = FALSE,
                parallel = TRUE,
                cores = n_cores)
print(mdri)

# let's look at naive FRR to get a sense
frr_untreated <- frrcal(data = untreated_all,
                        subid_var = "subject_id",
                        time_var = "days_since_eddi",
                        recency_cutoff_time = 730.5,
                        recency_rule = "independent_thresholds",
                        recency_vars = c("LAg_ODn","viral_load"),
                        recency_params = c(1.5,0,75,1),
                        debug = TRUE)
print(frr_untreated)
frr_treated <- frrcal(data = treated_all,
                      subid_var = "subject_id",
                      time_var = "days_since_eddi",
                      recency_cutoff_time = 730.5,
                      recency_rule = "independent_thresholds",
                      recency_vars = c("LAg_ODn","viral_load"),
                      recency_params = c(1.5,0,75,1))
print(frr_treated)