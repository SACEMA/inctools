#!/usr/bin/env Rscript
# Installation script for InctoolsJulia R package

cat("InctoolsJulia - R Interface to Inctools.jl\n")
cat("Package Installation\n")
cat(rep("=", 70), "\n\n", sep = "")

# Check R version
cat("Checking R version...\n")
r_version <- as.numeric(paste(R.version$major, R.version$minor, sep = "."))
if (r_version < 3.5) {
  stop("R version 3.5.0 or higher is required. You have: ", R.version.string)
}
cat("  OK: ", R.version.string, "\n\n")

# Install JuliaCall if needed
cat("Checking for JuliaCall package...\n")
if (!requireNamespace("JuliaCall", quietly = TRUE)) {
  cat("  Installing JuliaCall...\n")
  install.packages("JuliaCall")
  cat("  OK: JuliaCall installed\n\n")
} else {
  cat("  OK: JuliaCall already installed\n\n")
}

# Check Julia installation
cat("Checking for Julia...\n")
julia_path <- Sys.which("julia")
if (julia_path == "") {
  cat("  WARNING: Julia not found in PATH\n")
  cat("  Please install Julia from: https://julialang.org/downloads/\n")
  cat("  After installation, make sure 'julia' is in your PATH\n\n")
  cat("  To test if Julia is accessible, run in terminal:\n")
  cat("    julia --version\n\n")
  stop("Julia not found. Please install Julia and try again.")
} else {
  # Get Julia version
  julia_version <- system("julia --version", intern = TRUE)
  cat("  OK: ", julia_version, "\n")
  cat("  Path: ", julia_path, "\n\n")
}

# Initialize JuliaCall to verify it works
cat("Testing JuliaCall connection to Julia...\n")
tryCatch({
  julia <- JuliaCall::julia_setup()
  julia_ver <- JuliaCall::julia_eval("VERSION")
  cat("  OK: Successfully connected to Julia v", as.character(julia_ver), "\n\n")
}, error = function(e) {
  cat("  ERROR: Failed to connect to Julia\n")
  cat("  Error message: ", e$message, "\n\n")
  stop("JuliaCall setup failed")
})

# Check if Inctools.jl package exists
cat("Checking for Inctools.jl package...\n")
# Try to find Inctools in common locations
possible_paths <- c(
  file.path(getwd(), "..", "Inctools"),  # From InctoolsJulia/ directory
  file.path(getwd(), "Inctools"),        # From julia/ directory
  file.path(dirname(getwd()), "Inctools") # From subdirectory
)

inctools_path <- NULL
for (path in possible_paths) {
  if (dir.exists(path)) {
    inctools_path <- normalizePath(path)
    break
  }
}

if (is.null(inctools_path)) {
  cat("  ERROR: Inctools.jl directory not found\n")
  cat("  Tried locations:\n")
  for (p in possible_paths) {
    cat("    - ", p, "\n")
  }
  cat("\n  Please run this script from the julia/ or julia/InctoolsJulia/ directory\n\n")
  stop("Inctools.jl package not found")
}
cat("  OK: Found Inctools.jl at: ", inctools_path, "\n\n")

# Activate and precompile Inctools.jl
cat("Activating and precompiling Inctools.jl...\n")
cat("  (This may take 1-2 minutes on first run)\n")
tryCatch({
  JuliaCall::julia_eval(sprintf('using Pkg; Pkg.activate("%s")', inctools_path))
  JuliaCall::julia_eval("using Inctools")
  cat("  OK: Inctools.jl loaded successfully\n\n")
}, error = function(e) {
  cat("  ERROR: Failed to load Inctools.jl\n")
  cat("  Error message: ", e$message, "\n\n")
  stop("Inctools.jl loading failed")
})

cat(rep("=", 70), "\n", sep = "")
cat("INSTALLATION SUCCESSFUL!\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Next steps:\n")
cat("  1. To use the package, source the R files (from InctoolsJulia/ directory):\n")
cat("       source('R/zzz.R')\n")
cat("       source('R/inctools.R')\n\n")
cat("  2. Or install as an R package (from julia/ directory):\n")
cat("       install.packages('InctoolsJulia', repos = NULL, type = 'source')\n\n")
cat("  3. Run the test script (from julia/ directory):\n")
cat("       Rscript tests/test_R_api.R\n\n")
cat("  4. See InctoolsJulia/README.md for usage examples\n\n")

cat("Happy analyzing!\n")
