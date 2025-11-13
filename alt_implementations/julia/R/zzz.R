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
    pkg_path <- system.file("Inctools", package = "InctoolsJulia", mustWork = FALSE)
    if (pkg_path == "") {
      # If not installed as R package, try relative path
      pkg_path <- file.path(getwd(), "Inctools")
      if (!dir.exists(pkg_path)) {
        stop("Cannot find Inctools.jl package directory")
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
