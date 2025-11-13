# InctoolsJulia - R Interface to Inctools.jl

This package provides an R interface to the **Inctools.jl** Julia package for HIV incidence estimation using cross-sectional survey data with recency biomarkers.

**Note:** This package is named **InctoolsJulia** to distinguish it from the existing **inctools** package on CRAN. Both packages provide similar functionality, but this version uses a high-performance Julia backend.

## Features

- **Fast computation** using Julia backend with Gibbs sampling
- **Kassanjee method** for incidence estimation
- **Bootstrap and Delta method** uncertainty estimation
- **Simple R API** that feels native to R users
- **Automatic Julia setup** - just install and use

## Installation

### Prerequisites

1. **Install Julia** (version ≥ 1.6)
   - Download from: https://julialang.org/downloads/
   - Add Julia to your PATH

2. **Install JuliaCall R package**
   ```r
   install.packages("JuliaCall")
   ```

### Install InctoolsJulia for R

#### Option 1: Install from source (current directory)

```r
# From R, in the directory containing this README
install.packages(".", repos = NULL, type = "source")
```

#### Option 2: Load without installing

```r
# Load functions directly
source("R/zzz.R")
source("R/inctools.R")

# Initialize Julia
inctools_setup()
```

## Quick Start

```r
library(InctoolsJulia)

# Example 1: Estimate incidence from count data
result <- inccounts(
  n = 1000,           # Total sample size
  npos = 200,         # Number HIV positive
  ntestR = 180,       # Number tested for recency
  nR = 20,            # Number recent infections
  mdri = 130,         # Mean duration of recent infection (days)
  frr = 0.01,         # False recent rate
  sigma_mdri = 15,    # Standard error of MDRI
  sigma_frr = 0.005,  # Standard error of FRR
  bs = 2000           # Bootstrap samples
)

print(result)
# $I
# [1] 0.07525
#
# $CI
# [1] 0.03854 0.12345
#
# $RSE
# [1] 0.2842
```

## Usage Examples

### Example 1: Basic Incidence Estimation (Delta Method)

```r
library(InctoolsJulia)

# Using Delta method (faster, analytical)
result <- inccounts(
  n = 1000,
  npos = 200,
  ntestR = 180,
  nR = 20,
  mdri = 130,
  frr = 0.01,
  bs = 0  # 0 = Delta method
)

cat(sprintf("Incidence: %.4f (95%% CI: %.4f - %.4f)\n",
            result$I, result$CI[1], result$CI[2]))
```

### Example 2: Bootstrap Uncertainty Estimation

```r
# Using bootstrap (more robust for small samples)
result <- inccounts(
  n = 1000,
  npos = 200,
  ntestR = 180,
  nR = 20,
  mdri = 130,
  frr = 0.01,
  sigma_mdri = 15,
  sigma_frr = 0.005,
  bs = 2000  # 2000 bootstrap samples
)
```

### Example 3: Accounting for Covariance

```r
# When prevalence and recent prevalence are correlated
result <- inccounts(
  n = 1000,
  npos = 200,
  ntestR = 180,
  nR = 20,
  mdri = 130,
  frr = 0.01,
  sigma_mdri = 15,
  sigma_frr = 0.005,
  covar = 0.0002,  # Covariance between prev and prevR
  bs = 2000
)
```

### Example 4: Using Prevalence Proportions Directly

```r
# If you already have prevalence estimates
result <- incprops(
  prev = 0.20,          # Prevalence
  sigma_prev = 0.015,   # SE of prevalence
  prevR = 0.11,         # Recent infection prevalence
  sigma_prevR = 0.023,  # SE of recent prevalence
  mdri = 130,
  sigma_mdri = 15,
  frr = 0.01,
  sigma_frr = 0.005,
  bs = 2000
)
```

### Example 5: Comparing Two Populations

```r
# Test if incidence differs between two groups
result <- incdif(
  prev = c(0.20, 0.15),
  sigma_prev = c(0.015, 0.012),
  prevR = c(0.10, 0.08),
  sigma_prevR = c(0.02, 0.015),
  mdri = 130,
  sigma_mdri = 15,
  frr = 0.01,
  sigma_frr = 0.005,
  bs = 2000
)

cat(sprintf("Difference: %.4f (95%% CI: %.4f - %.4f), p = %.4f\n",
            result$Delta, result$CI[1], result$CI[2], result$p))
```

### Example 6: Simple Prevalence Estimation

```r
# Calculate prevalence with confidence interval
result <- prevalence(
  pos = 100,
  n = 1000,
  ci = TRUE
)

print(result)
```

## Function Reference

### Core Functions

- **`inccounts()`** - Estimate incidence from count data
- **`incprops()`** - Estimate incidence from prevalence proportions
- **`incdif()`** - Test difference in incidence between populations
- **`prevalence()`** - Calculate prevalence with confidence intervals

### Utility Functions

- **`inctools_setup()`** - Manually initialize Julia (usually automatic)

## Parameters

### Common Parameters

- **`mdri`** - Mean duration of recent infection (in days)
- **`frr`** - False recent rate (proportion, 0-1)
- **`T`** - Time cutoff (default: 730.5 days = 2 years)
- **`bs`** - Number of bootstrap samples (0 = Delta method)
- **`alpha`** - Significance level (default: 0.05)
- **`covar`** - Covariance between prevalence estimates (default: 0)

### Design Effects

- **`de_npos`** - Design effect for prevalence (default: 1.0)
- **`de_nR`** - Design effect for recent prevalence (default: 1.0)

## Performance Notes

The Julia backend uses:
- **Gibbs sampling** for independent variables (very fast)
- **Rejection sampling** when covariance > 0 (slower, handles correlation)

First call will be slow (~30 seconds) due to Julia compilation. Subsequent calls are fast (<1 second).

## Troubleshooting

### Julia not found

```r
# Specify Julia path explicitly
inctools_setup(julia_path = "/path/to/julia/bin")
```

### Package loading errors

```r
# Check Julia installation
JuliaCall::julia_setup()

# Check Julia can find Inctools.jl
JuliaCall::julia_eval('using Pkg; Pkg.activate("./Inctools"); using Inctools')
```

### Slow first run

This is normal! Julia compiles code on first use. Subsequent runs are much faster.

## Comparison with R inctools Package

| Feature | Inctools.jl (this) | inctools (R) |
|---------|-------------------|--------------|
| Language | Julia (via R) | Pure R |
| Speed | Very fast | Moderate |
| Bootstrap | ✓ | ✓ |
| Delta method | ✓ | ✓ |
| Gibbs sampling | ✓ | ✗ |
| Covariance support | ✓ | ✓ |
| Installation | Requires Julia | R only |

## Citation

If you use this package, please cite:

```
Kassanjee R, et al. (2012). A new general biomarker-based incidence estimator.
Epidemiology, 23(5), 721-728.
```

## License

MIT License

## Authors

- Eduard Grebe (original Julia implementation)
- Claude AI (R wrapper)

## Links

- GitHub: https://github.com/SACEMA/inctools
- Julia Package: See `Inctools/` directory
- R inctools: https://cran.r-project.org/package=inctools
