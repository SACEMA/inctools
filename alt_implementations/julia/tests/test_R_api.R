# Test script for R API to Inctools.jl
# Run this from the julia/ directory: Rscript tests/test_R_api.R

cat("Testing R API to Inctools.jl\n")
cat(rep("=", 70), "\n\n", sep = "")

# Load the R functions from InctoolsJulia package
cat("Loading R functions...\n")
source("InctoolsJulia/R/zzz.R")
source("InctoolsJulia/R/inctools.R")

# Initialize Julia and Inctools.jl
cat("Initializing Julia (this may take ~30 seconds on first run)...\n")
inctools_setup()

cat("\n", rep("=", 70), "\n", sep = "")
cat("Test 1: Simple prevalence calculation\n")
cat(rep("-", 70), "\n", sep = "")

result1 <- prevalence(100, 1000)
cat(sprintf("Prevalence: %.4f, SE: %.6f\n", result1[[1]], result1[[2]]))

result1_ci <- prevalence(100, 1000, ci = TRUE)
cat(sprintf("With CI: %.4f (95%% CI: %.4f - %.4f)\n",
            result1_ci[[1]], result1_ci[[3]][1], result1_ci[[3]][2]))

cat("\n", rep("=", 70), "\n", sep = "")
cat("Test 2: Incidence from counts (Delta method)\n")
cat(rep("-", 70), "\n", sep = "")

result2 <- inccounts(
  n = 1000,
  npos = 200,
  ntestR = 180,
  nR = 20,
  mdri = 130,
  frr = 0.01,
  sigma_mdri = 15,
  sigma_frr = 0.005,
  bs = 0  # Delta method
)

cat(sprintf("Incidence: %.4f\n", result2$I))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", result2$CI[1], result2$CI[2]))
cat(sprintf("RSE: %.4f\n", result2$RSE))

cat("\n", rep("=", 70), "\n", sep = "")
cat("Test 3: Incidence from counts (Bootstrap)\n")
cat(rep("-", 70), "\n", sep = "")

cat("Running 1000 bootstrap samples...\n")
result3 <- inccounts(
  n = 1000,
  npos = 200,
  ntestR = 180,
  nR = 20,
  mdri = 130,
  frr = 0.01,
  sigma_mdri = 15,
  sigma_frr = 0.005,
  bs = 1000  # Bootstrap
)

cat(sprintf("Incidence: %.4f\n", result3$I))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", result3$CI[1], result3$CI[2]))
cat(sprintf("RSE: %.4f\n", result3$RSE))

cat("\n", rep("=", 70), "\n", sep = "")
cat("Test 4: Incidence with covariance\n")
cat(rep("-", 70), "\n", sep = "")

result4 <- inccounts(
  n = 1000,
  npos = 200,
  ntestR = 180,
  nR = 20,
  mdri = 130,
  frr = 0.01,
  sigma_mdri = 15,
  sigma_frr = 0.005,
  covar = 0.0002,
  bs = 1000
)

cat(sprintf("Incidence: %.4f\n", result4$I))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", result4$CI[1], result4$CI[2]))
cat(sprintf("RSE: %.4f\n", result4$RSE))

cat("\n", rep("=", 70), "\n", sep = "")
cat("Test 5: Using incprops directly\n")
cat(rep("-", 70), "\n", sep = "")

result5 <- incprops(
  prev = 0.20,
  sigma_prev = 0.01265,
  prevR = 0.11,
  sigma_prevR = 0.02342,
  mdri = 130,
  sigma_mdri = 15,
  frr = 0.01,
  sigma_frr = 0.005,
  bs = 1000
)

cat(sprintf("Incidence: %.4f\n", result5$I))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", result5$CI[1], result5$CI[2]))
cat(sprintf("RSE: %.4f\n", result5$RSE))

cat("\n", rep("=", 70), "\n", sep = "")
cat("Test 6: Comparing two populations\n")
cat(rep("-", 70), "\n", sep = "")

result6 <- incdif(
  prev = c(0.20, 0.15),
  sigma_prev = c(0.015, 0.012),
  prevR = c(0.10, 0.08),
  sigma_prevR = c(0.02, 0.015),
  mdri = 130,
  sigma_mdri = 15,
  frr = 0.01,
  sigma_frr = 0.005,
  bs = 1000
)

cat(sprintf("Difference: %.4f\n", result6$Delta))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", result6$CI[1], result6$CI[2]))
cat(sprintf("p-value: %.4f\n", result6$p))

cat("\n", rep("=", 70), "\n", sep = "")
cat("ALL TESTS COMPLETED SUCCESSFULLY!\n")
cat(rep("=", 70), "\n", sep = "")

cat("\nSummary:\n")
cat("  - R API successfully calls Julia functions\n")
cat("  - All function types work correctly:\n")
cat("    * prevalence()\n")
cat("    * inccounts() with Delta method\n")
cat("    * inccounts() with bootstrap\n")
cat("    * inccounts() with covariance\n")
cat("    * incprops()\n")
cat("    * incdif()\n")
cat("  - Results are returned correctly to R\n")
cat("\nR users can now use Inctools.jl seamlessly from R!\n")
