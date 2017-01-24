
ss_baseline_cohort(Inc=0.02, # Incidence at start
                   Prev=0.4, # Prevalence at start
                   FracIncRed=0.5, # Incidence reduction as proportion
                   n = "out", # ask for sample size to achieve specified power
                   Power = 0.8, # desired power to detect incidence reduction
                   alpha = 0.05, # confidence level
                   MDRI=200, 
                   RSE_MDRI=0.05, 
                   FRR=0.005, 
                   RSE_FRR=0.25, 
                   CR = 1, # coverage rate of recency test in baseline survey
                   DE_H = 1, # design effect on HIV prevalence in baseline survey
                   DE_R = 1, # design effect on prevalence of recency in baseline survey
                   DE_C = 1, # design effect on cohort incidence
                   BigT = 730, # T in days
                   CohortCR = 1, # proportion of negatives recruited into cohort
                   FUT = 1) # cohort follow-up time in same units as incidence (usually years)


prepimpact(Inc=0.01, 
           Prev=0.2, 
           FracIncRed=0.5, 
           n = "out", 
           Power = 0.8, 
           alpha = 0.05, 
           MDRI=180, 
           RSE_MDRI=0.1, 
           FRR=0.005, 
           RSE_FRR=0.25, 
           CR = 1, 
           DE_H = 2,
           DE_R = 2, 
           DE_C = 2, 
           BigT = 730, 
           CohortCR = 1, 
           FUT = 1.5)

# Example where desired power cannot be achieved
prepimpact(Inc=0.005, 
           Prev=0.1, 
           FracIncRed=0.25, 
           n = "out", 
           Power = 0.8, 
           alpha = 0.05, 
           MDRI=180, 
           RSE_MDRI=0.1, 
           FRR=0.005, 
           RSE_FRR=0.25, 
           CR = 1, 
           DE_H = 1.5,
           DE_R = 1.5, 
           DE_C = 1.5, 
           BigT = 730, 
           CohortCR = 1, 
           FUT = 1)
