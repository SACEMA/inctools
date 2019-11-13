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
        cores = 4,
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
        cores = 4,
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

inccounts(N = 10000,
          N_H = 1000,
          N_testR = 990,
          N_R = 99,
          DE_H = 1.2,
          DE_R = 1.2,
          Boot = FALSE,
          alpha = 0.05,
          MDRI= 200,
          RSE_MDRI = 0.1,
          FRR = 0.01,
          RSE_FRR = 0.25,
          BigT = 730.5,
          Covar_HR = 0,
          debug = FALSE)

incprops(PrevH = 0.1, 
         RSE_PrevH = 0.03286335,
         PrevR = 0.1, 
         RSE_PrevR = 0.1044466,
         Boot = FALSE,
         alpha = 0.05,
         MDRI= 200,
         RSE_MDRI = 0.1,
         FRR = 0.01,
         RSE_FRR = 0.25,
         BigT = 730.5,
         Covar_HR = 0,
         debug = FALSE)


incpower(I1 = 0.05, 
         I2 = 0.03, 
         PrevH1 = 0.20, 
         PrevH2 = 0.20,
         n1 = 5000, 
         n2 = 5000, 
         alpha = 0.05, 
         Power = "out", 
         SS = NULL,
         DE_H = c(1,1.1), 
         DE_R = 1, 
         BMest = "FRR.indep", 
         MDRI = 200,
         RSE_MDRI = 0.05, 
         FRR = c(0.01, 0.02), 
         RSE_FRR = c(0.20, 0.15), 
         BigT = 730.5)$Inc.Difference.Statistics$Power


