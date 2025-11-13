# Package initialization
.onLoad <- function(libname, pkgname) {
  # Check if JuliaCall is available
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("Package 'JuliaCall' is required but not installed. Install it with: install.packages('JuliaCall')")
  }

  packageStartupMessage("Loading Inctools.jl Julia package...")
  packageStartupMessage("Note: First run may take time to precompile Julia packages.")
}

# Initialize Julia and load Inctools package
.inctools_julia_setup <- function() {
  if (!exists(".inctools_julia_initialized", envir = .GlobalEnv)) {
    # Initialize Julia
    julia <- JuliaCall::julia_setup()

    # Get the path to the Inctools.jl package
    # When installed: inst/Inctools becomes <package-root>/Inctools
    pkg_path <- system.file("Inctools", package = "InctoolsJulia", mustWork = FALSE)

    if (pkg_path == "") {
      # If not installed as R package, try common development locations
      # 1. Try ../Inctools (from InctoolsJulia/ directory)
      # 2. Try ../../Inctools (from InctoolsJulia/R/ when sourcing)
      possible_paths <- c(
        file.path(getwd(), "..", "Inctools"),
        file.path(getwd(), "Inctools"),
        file.path(dirname(getwd()), "Inctools")
      )

      for (path in possible_paths) {
        if (dir.exists(path)) {
          pkg_path <- normalizePath(path)
          break
        }
      }

      if (pkg_path == "" || !dir.exists(pkg_path)) {
        stop("Cannot find Inctools.jl package directory. ",
             "Make sure you are in the julia/ directory or install the package properly.")
      }
    }

    # Activate the Inctools.jl project
    JuliaCall::julia_eval(sprintf('using Pkg; Pkg.activate("%s")', pkg_path))

    # Load Inctools.jl
    JuliaCall::julia_eval("using Inctools")

    assign(".inctools_julia_initialized", TRUE, envir = .GlobalEnv)
    message("Inctools.jl loaded successfully!")
  }
}
