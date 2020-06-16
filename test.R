# Check build using R-hub
# Interactive
rhub::check(path = "~/dev/inctools/inctools/")
# Automated for CRAN
rhub::check_for_cran(path = "~/dev/inctools/inctools/")

rhub::check_with_roldrel(path = "~/dev/inctools/inctools/")
rhub::check_with_rrelease(path = "~/dev/inctools/inctools/")
rhub::check_with_rpatched(path = "~/dev/inctools/inctools/")
rhub::check_with_rdevel(path = "~/dev/inctools/inctools/") 


n_cores <- 12
setwd("~/dev/incidence-kzn/")

library(dplyr); library(tidyr); library(inctools);

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

# Test whether bsparms are output

mdrical(data = untreated_c,
        subid_var = "subject_id",
        time_var = "days_since_eddi",
        functional_forms = "logit_cubic",
        recency_cutoff_time = 730.5,
        recency_rule = "independent_thresholds",
        recency_vars = c("LAg_ODn","viral_load"),
        recency_params = c(1.5,0,75,1),
        n_bootstraps = 10,
        plot = TRUE,
        parallel = TRUE,
        cores = 6,
        output_bs_parms = TRUE,
        debug = FALSE)

mdrical(data = untreated_c,
        subid_var = "subject_id",
        time_var = "days_since_eddi",
        functional_forms = "cloglog_linear",
        recency_cutoff_time = 730.5,
        recency_rule = "independent_thresholds",
        recency_vars = c("LAg_ODn","viral_load"),
        recency_params = c(1.5,0,75,1),
        n_bootstraps = 10,
        plot = TRUE,
        parallel = TRUE,
        cores = 6,
        output_bs_parms = TRUE,
        debug = FALSE)

mdrical(data = untreated_c,
        subid_var = "subject_id",
        time_var = "days_since_eddi",
        functional_forms = "logit_cubic",
        recency_cutoff_time = 730.5,
        recency_rule = "independent_thresholds",
        recency_vars = c("LAg_ODn","viral_load"),
        recency_params = c(1.5,0,75,1),
        n_bootstraps = 10,
        plot = TRUE,
        parallel = FALSE,
        output_bs_parms = TRUE,
        debug = FALSE)

mdrical(data = untreated_c,
        subid_var = "subject_id",
        time_var = "days_since_eddi",
        functional_forms = "cloglog_linear",
        recency_cutoff_time = 730.5,
        recency_rule = "independent_thresholds",
        recency_vars = c("LAg_ODn","viral_load"),
        recency_params = c(1.5,0,75,1),
        n_bootstraps = 10,
        plot = TRUE,
        parallel = FALSE,
        output_bs_parms = TRUE,
        debug = FALSE)

mdrical(data = untreated_c,
        subid_var = "subject_id",
        time_var = "days_since_eddi",
        functional_forms = "cloglog_linear",
        recency_cutoff_time = 730.5,
        recency_rule = "independent_thresholds",
        recency_vars = c("LAg_ODn","viral_load"),
        recency_params = c(1.5,0,75,1),
        n_bootstraps = 0,
        plot = TRUE,
        debug = FALSE)

# with both functional forms
mdrical(data = untreated_c,
        subid_var = "subject_id",
        time_var = "days_since_eddi",
        functional_forms = c("cloglog_linear","logit_cubic"),
        recency_cutoff_time = 730.5,
        recency_rule = "independent_thresholds",
        recency_vars = c("LAg_ODn","viral_load"),
        recency_params = c(1.5,0,75,1),
        n_bootstraps = 10,
        plot = TRUE,
        parallel = FALSE,
        output_bs_parms = TRUE,
        debug = FALSE)

mdrical(data = untreated_c,
        subid_var = "subject_id",
        time_var = "days_since_eddi",
        functional_forms = c("cloglog_linear","logit_cubic"),
        recency_cutoff_time = 730.5,
        recency_rule = "independent_thresholds",
        recency_vars = c("LAg_ODn","viral_load"),
        recency_params = c(1.5,0,75,1),
        n_bootstraps = 10,
        plot = TRUE,
        parallel = TRUE,
        cores = 4,
        output_bs_parms = TRUE,
        debug = FALSE)



# Other tests


mdri <- mdrical(data = untreated_c,
                subid_var = "subject_id",
                time_var = "days_since_eddi",
                functional_forms = "logit_cubic",
                recency_cutoff_time = 730.5,
                recency_rule = "independent_thresholds",
                recency_vars = c("LAg_ODn","viral_load"),
                recency_params = c(1.5,0,75,1),
                n_bootstraps = 100,
                plot = TRUE,
                parallel = TRUE,
                output_bs_parms = TRUE,
                cores = n_cores,
                debug = TRUE)
print(mdri)

# let's look at naive FRR to get a sense
frr_untreated <- frrcal(data = untreated_all,
                        subid_var = "subject_id",
                        time_var = "days_since_eddi",
                        recency_cutoff_time = 730.5,
                        recency_rule = "independent_thresholds",
                        recency_vars = c("LAg_ODn","viral_load"),
                        recency_params = c(1.5,0,75,1),
                        method = "exact",
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

I1 <- incprops(PrevH = 0.20, RSE_PrevH = 0.028, PrevR = 0.10, RSE_PrevR = 0.09,
         BS_Count = 10000, Boot = TRUE, MDRI = 200, RSE_MDRI = 0.05,
         FRR = 0.01, RSE_FRR = 0.2, BigT = 730, debug = TRUE)
print(I1)
str(I1)

I2 <- incprops(PrevH = c(0.20,0.21,0.18), 
               RSE_PrevH = c(0.028,0.03,0.022),
               PrevR = c(0.10,0.13,0.12), 
               RSE_PrevR = c(0.094,0.095,0.05),
               Boot = FALSE, 
               BMest = 'MDRI.FRR.indep',
               MDRI = c(200,180,180), 
               RSE_MDRI = c(0.05,0.07,0.06),
               FRR = c(0.01,0.009,0.02), 
               RSE_FRR = c(0.2,0.2,0.1), 
               BigT = 730.5)
I2$Incidence.Difference.Statistics$p.value
print(I2)
str(I2)

c(mdri$MDRI$SE, mdri$MDRI$CI_LB, mdri$MDRI$CI_UB, 
  nrow(mdri$BSparms$logit_cubic),
  mean(mdri$BSparms$logit_cubic$beta0),
  mean(mdri$BSparms$logit_cubic$beta1))




SEs <- 0
LBs <- 0
UBs <- 0
for (i in 1:25) {
  mdri <-mdrical(data=excalibdata,
                    subid_var = "SubjectID",
                    time_var = "DaysSinceEDDI",
                    recency_cutoff_time = 730.5,
                    inclusion_time_threshold = 800,
                    functional_forms = c("cloglog_linear"),
                    recency_rule = "binary_data",
                    recency_vars = "Recent",
                    n_bootstraps = 250,
                    random_seed = 123,
                    plot = FALSE,
                    parallel = TRUE)
  SEs[i] <- mdri$MDRI$SE
  LBs[i] <- mdri$MDRI$CI_LB
  UBs[i] <- mdri$MDRI$CI_UB
}
min(SEs) - 0.5
max(LBs) + 0.5
min(UBs) - 0.5

SEs <- 0
LBs <- 0
UBs <- 0
for (i in 1:25) {
  mdri <- mdrical(data=excalibdata,
                  subid_var = "SubjectID",
                  time_var = "DaysSinceEDDI",
                  recency_cutoff_time = 730.5,
                  inclusion_time_threshold = 800,
                  functional_forms = c("logit_cubic"),
                  recency_rule = "binary_data",
                  recency_vars = "Recent",
                  n_bootstraps = 250,
                  random_seed = 123,
                  plot = FALSE,
                  parallel = TRUE,
                  cores = 2)
  SEs[i] <- mdri$MDRI$SE
  LBs[i] <- mdri$MDRI$CI_LB
  UBs[i] <- mdri$MDRI$CI_UB
}
min(SEs) - 0.5
max(LBs) + 0.5
min(UBs) - 0.5


