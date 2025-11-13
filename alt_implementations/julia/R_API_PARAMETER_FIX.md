# R API Parameter Name Fix

**Date:** 2025-11-13
**Issue:** R wrapper was using Latin parameter names (alpha, sigma_mdri, sigma_frr) but Julia functions expect Greek letters (α, σ_mdri, σ_frr)
**Status:** ✅ Fixed

---

## Problem

When calling Julia functions from R via JuliaCall, the **keyword parameter names** must exactly match the Julia function signatures. The Julia Inctools package uses Greek letters for some parameters:

- `α` (Greek alpha) instead of `alpha`
- `σ_mdri` (Greek sigma) instead of `sigma_mdri`
- `σ_frr` (Greek sigma) instead of `sigma_frr`

### Error Example

```r
> result1 <- prevalence(100, 1000)
Error in Julia:
MethodError: no method matching prevalence(::Float64, ::Float64, ::Float64;
                                           ci::Bool, alpha::Float64)
This method does not support all of the given keyword arguments.
Closest candidates are:
  prevalence(::Any, ::Any, ::Any; ci, α) got unsupported keyword argument "alpha"
```

---

## Solution

Updated all R wrapper functions in `R/inctools.R` to use Greek letters when passing keyword arguments to Julia, while keeping the R function signatures user-friendly with Latin letters.

### Pattern

```r
# R function signature (user-friendly)
my_function <- function(..., alpha = 0.05, sigma_mdri = 0.0) {

  # Call Julia with Greek letters for keyword arguments
  result <- JuliaCall::julia_call("my_function",
                                   ...,
                                   α = alpha,        # Convert here
                                   σ_mdri = sigma_mdri)  # Convert here
}
```

---

## Changes Made

### 1. `prevalence()` function

**Before:**
```r
result <- JuliaCall::julia_call("prevalence", pos, n, de,
                                 ci = ci, alpha = alpha)
```

**After:**
```r
result <- JuliaCall::julia_call("prevalence", pos, n, de,
                                 ci = ci, α = alpha)  # ✓ Changed
```

### 2. `incprops()` function

**Before:**
```r
result <- JuliaCall::julia_call("incprops",
                                 prev, sigma_prev, prevR, sigma_prevR,
                                 mdri, sigma_mdri, frr, sigma_frr,
                                 covar = covar, T = T,
                                 timeconversion = timeconversion,
                                 bs = bs, alpha = alpha, per = per)
```

**After:**
```r
result <- JuliaCall::julia_call("incprops",
                                 prev, sigma_prev, prevR, sigma_prevR,
                                 mdri, sigma_mdri, frr, sigma_frr,
                                 covar = covar, T = T,
                                 timeconversion = timeconversion,
                                 bs = bs, α = alpha, per = per)  # ✓ Changed
```

**Note:** The positional parameters `sigma_prev`, `sigma_prevR`, etc. don't need to match the Julia names because they're passed by position, not by name.

### 3. `inccounts()` function

**Before:**
```r
result <- JuliaCall::julia_call("inccounts",
                                 n, npos, ntestR, nR, mdri, frr,
                                 de_npos = de_npos, de_nR = de_nR,
                                 sigma_mdri = sigma_mdri, sigma_frr = sigma_frr,
                                 covar = covar, T = T,
                                 timeconversion = timeconversion,
                                 bs = bs, alpha = alpha, per = per)
```

**After:**
```r
result <- JuliaCall::julia_call("inccounts",
                                 n, npos, ntestR, nR, mdri, frr,
                                 de_npos = de_npos, de_nR = de_nR,
                                 σ_mdri = sigma_mdri, σ_frr = sigma_frr,  # ✓ Changed
                                 covar = covar, T = T,
                                 timeconversion = timeconversion,
                                 bs = bs, α = alpha, per = per)  # ✓ Changed
```

### 4. `incdif()` function

**Before:**
```r
result <- JuliaCall::julia_call("incdif",
                                 prev, sigma_prev, prevR, sigma_prevR,
                                 mdri, sigma_mdri, frr, sigma_frr,
                                 covar = covar, T = T,
                                 timeconversion = timeconversion,
                                 bs = bs, alpha = alpha, per = per)
```

**After:**
```r
result <- JuliaCall::julia_call("incdif",
                                 prev, sigma_prev, prevR, sigma_prevR,
                                 mdri, sigma_mdri, frr, sigma_frr,
                                 covar = covar, T = T,
                                 timeconversion = timeconversion,
                                 bs = bs, α = alpha, per = per)  # ✓ Changed
```

**Note:** For `incdif()`, the `sigma_*` parameters are positional (passed before the semicolon in Julia), so they don't need Greek letters. Only the keyword parameter `α` needs the Greek letter.

---

## Julia Function Signatures Reference

### `prevalence`
```julia
function prevalence(pos, n, de = 1; ci = false, α = 0.05)
```
- Keyword: `α` ✓

### `incprops`
```julia
function incprops(prev::Float64, σ_prev::Float64, prevR::Float64, σ_prevR::Float64,
                  mdri::Float64, σ_mdri::Float64, frr::Float64, σ_frr::Float64;
                  covar::Float64 = 0.0, T = 730.5, timeconversion = 365.25,
                  bs::Int64 = 0, gibbs = false, bs_numbers = false,
                  bs_numbers_n::AbstractVector{Int64} = [0, 0],
                  α::Float64 = 0.05, per::Int64 = 1)
```
- Positional: `σ_prev`, `σ_prevR`, `σ_mdri`, `σ_frr` (passed by position, names don't matter)
- Keyword: `α` ✓

### `inccounts`
```julia
function inccounts(n::Int64, npos::Int64, ntestR::Int64, nR::Int64,
                   mdri::Float64, frr::Float64;
                   de_npos::Float64 = 1.0, de_nR::Float64 = 1.0,
                   σ_mdri::Float64 = 0.0, σ_frr::Float64 = 0.0,
                   covar::Union{Float64, Nothing} = nothing,
                   T = 730.5, timeconversion = 365.25,
                   bs::Int64 = 0, α::Float64 = 0.05, per::Int64 = 1)
```
- Keywords: `σ_mdri`, `σ_frr`, `α` ✓

### `incdif`
```julia
function incdif(prev::AbstractVector{Float64}, σ_prev::AbstractVector{Float64},
                prevR::AbstractVector{Float64}, σ_prevR::AbstractVector{Float64},
                mdri::Float64, σ_mdri::Float64, frr::Float64, σ_frr::Float64;
                covar::AbstractVector{Float64}, T = 730.5, timeconversion = 365.25,
                bs::Int64 = 0, output_bs = false, bs_numbers = false,
                bs_numbers_n::AbstractVector{Int64} = [0, 0, 0, 0],
                α::Float64 = 0.05, bonf_cor::Int64 = 1, per::Int64 = 1)
```
- Positional: `σ_prev`, `σ_prevR`, `σ_mdri`, `σ_frr` (passed by position, names don't matter)
- Keyword: `α` ✓

---

## Key Distinction: Positional vs Keyword Parameters

### Positional Parameters (before semicolon in Julia)
- Names **don't matter** in JuliaCall
- Passed by **order**, not by name
- Example: `sigma_prev` in R → `σ_prev` in Julia (automatic, no conversion needed)

### Keyword Parameters (after semicolon in Julia)
- Names **must match exactly**
- Must use Greek letters if Julia uses them
- Example: `alpha = 0.05` in R → `α = 0.05` when calling Julia

---

## User Impact

### For R Users (NO CHANGE)

R users continue to use familiar Latin letters:

```r
# R code looks exactly the same
result <- prevalence(100, 1000, alpha = 0.05)
result <- inccounts(1000, 200, 180, 20, 130, 0.01,
                    sigma_mdri = 15, sigma_frr = 0.005, alpha = 0.05)
```

### Behind the Scenes

The R wrapper automatically converts to Greek letters when calling Julia:

```r
prevalence <- function(pos, n, de = 1.0, ci = FALSE, alpha = 0.05) {
  # User passes alpha = 0.05
  result <- JuliaCall::julia_call("prevalence", pos, n, de,
                                   ci = ci,
                                   α = alpha)  # Automatically converted to α
  return(result)
}
```

---

## Testing

All functions now work correctly:

```r
source("R/zzz.R")
source("R/inctools.R")
inctools_setup()

# All of these now work
prevalence(100, 1000)
prevalence(100, 1000, ci = TRUE, alpha = 0.05)
incprops(0.20, 0.015, 0.10, 0.02, 130, 15, 0.01, 0.005, alpha = 0.05)
inccounts(1000, 200, 180, 20, 130, 0.01,
          sigma_mdri = 15, sigma_frr = 0.005, alpha = 0.05)
incdif(c(0.20, 0.15), c(0.015, 0.012), c(0.10, 0.08), c(0.02, 0.015),
       130, 15, 0.01, 0.005, alpha = 0.05)
```

---

## Summary

✅ **All R wrapper functions fixed**

The fix was simple but critical:
- R users use `alpha`, `sigma_mdri`, `sigma_frr` (easy to type)
- R wrapper converts to `α`, `σ_mdri`, `σ_frr` when calling Julia
- No change needed to user-facing API
- All functions now work correctly

**Files Modified:**
- `R/inctools.R` (4 functions updated)

**Status:** Ready for use
