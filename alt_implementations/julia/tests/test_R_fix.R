# Quick test of the R API fix
# Run from julia/ directory: Rscript tests/test_R_fix.R
source("InctoolsJulia/R/zzz.R")
source("InctoolsJulia/R/inctools.R")

cat("Initializing Julia...\n")
inctools_setup()

cat("\nTest 1: prevalence\n")
result1 <- prevalence(100, 1000)
print(result1)

cat("\nTest 2: prevalence with CI\n")
result2 <- prevalence(100, 1000, ci = TRUE)
print(result2)

cat("\n✓ Tests passed!\n")
