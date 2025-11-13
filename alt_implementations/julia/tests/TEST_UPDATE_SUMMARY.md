# Test Files Update for New Directory Structure

**Date:** 2025-11-13
**Status:** ✅ Complete

---

## Changes Made

Updated all test scripts to work with the new directory structure where:
- Julia package is in `Inctools/` subdirectory
- R package is in `InctoolsJulia/` subdirectory
- Test files are in `tests/` subdirectory
- Tests are run from `julia/` directory (or `tests/` for Julia)

---

## Files Updated

### Julia Test Files (4 files)

Added smart path detection to all Julia tests:

**Files:**
- `tests/test_comprehensive.jl`
- `tests/test_rtmvnorm.jl`
- `tests/test_gibbs_vs_rejection.jl`
- `tests/test_covar_deprecation.jl`

**Before:**
```julia
using Pkg
Pkg.activate("./Inctools")
using Inctools
```

**After:**
```julia
using Pkg

# Find Inctools package (works from julia/ or tests/ directory)
inctools_path = isdir("./Inctools") ? "./Inctools" : "../Inctools"
if !isdir(inctools_path)
    error("Cannot find Inctools package. Run this from julia/ directory.")
end
Pkg.activate(inctools_path)

using Inctools
```

**Result:** Julia tests now work from both `julia/` and `tests/` directories!

### R Test Files (2 files)

### 1. `tests/test_R_api.R`

**Before:**
```r
source("R/zzz.R")
source("R/inctools.R")
```

**After:**
```r
source("InctoolsJulia/R/zzz.R")
source("InctoolsJulia/R/inctools.R")
```

### 2. `tests/test_R_fix.R`

**Before:**
```r
source("R/zzz.R")
source("R/inctools.R")
```

**After:**
```r
source("InctoolsJulia/R/zzz.R")
source("InctoolsJulia/R/inctools.R")
```

---

## Usage

### Running Julia Tests

**Option 1: From julia/ directory**
```bash
julia tests/test_comprehensive.jl
julia tests/test_rtmvnorm.jl
julia tests/test_gibbs_vs_rejection.jl
julia tests/test_covar_deprecation.jl
```

**Option 2: From tests/ directory**
```bash
cd tests
julia test_comprehensive.jl
julia test_rtmvnorm.jl
julia test_gibbs_vs_rejection.jl
julia test_covar_deprecation.jl
```

### Running R Tests

**From julia/ directory (required):**
```bash
Rscript tests/test_R_api.R
Rscript tests/test_R_fix.R
```

---

## Directory Structure Expected by Tests

```
julia/                         # ← Run all tests from here
├── Inctools/                  # Julia package
│   ├── Project.toml
│   └── src/Inctools.jl
│
├── InctoolsJulia/             # R package
│   ├── R/
│   │   ├── zzz.R              # ← R tests source this
│   │   └── inctools.R         # ← R tests source this
│   ├── inst/
│   │   └── Inctools -> ../../Inctools
│   └── DESCRIPTION
│
└── tests/                     # All test files
    ├── test_comprehensive.jl
    ├── test_R_api.R           # ← Updated
    ├── test_R_fix.R           # ← Updated
    └── ...
```

---

## Path Resolution

### Julia Tests

**From julia/ directory:**
- Checks `./Inctools` → Found ✓
- Activates `./Inctools`

**From tests/ directory:**
- Checks `./Inctools` → Not found
- Checks `../Inctools` → Found ✓
- Activates `../Inctools`

### R Tests

**From julia/ directory (required):**
- `InctoolsJulia/R/zzz.R` → Resolves to correct file ✓
- `InctoolsJulia/R/inctools.R` → Resolves to correct file ✓

The R functions in `zzz.R` then find the Julia package:
- Checks `../Inctools` (from `InctoolsJulia/`)
- Finds `Inctools/Project.toml` ✓

---

## Verification

### Verify Julia Tests

```bash
# From julia/ directory
$ julia tests/test_comprehensive.jl
  Activating project at `.../julia/Inctools`
Comprehensive Inctools Testing
======================================================================
[9 tests pass]

# From tests/ directory
$ cd tests
$ julia test_comprehensive.jl
  Activating project at `.../julia/Inctools`
Comprehensive Inctools Testing
======================================================================
[9 tests pass]
```

### Verify R Tests

```bash
# From julia/ directory
$ Rscript tests/test_R_api.R
Testing R API to Inctools.jl
======================================================================
Loading R functions...
Initializing Julia...
[6 tests pass]
```

---

## Notes

- **Julia tests** now have smart path detection
  - Work from `julia/` directory (uses `./Inctools`)
  - Work from `tests/` directory (uses `../Inctools`)
  - Automatic detection with clear error if run from wrong location
- **R tests** must be run from `julia/` directory
  - They use relative paths to `InctoolsJulia/R/`
  - R's `zzz.R` also has smart path detection for finding Inctools

---

## Checklist

### Julia Tests
- [x] Added smart path detection to `test_comprehensive.jl`
- [x] Added smart path detection to `test_rtmvnorm.jl`
- [x] Added smart path detection to `test_gibbs_vs_rejection.jl`
- [x] Added smart path detection to `test_covar_deprecation.jl`
- [x] Verified tests work from `julia/` directory
- [x] Verified tests work from `tests/` directory

### R Tests
- [x] Updated `tests/test_R_api.R` source paths
- [x] Updated `tests/test_R_fix.R` source paths
- [x] Added usage comments in test files
- [x] Verified tests work from `julia/` directory

### Documentation
- [x] Created `tests/TEST_README.md`
- [x] Updated `tests/TEST_UPDATE_SUMMARY.md` (this file)
- [x] Documented all changes and usage patterns

---

## Summary

✅ **All test files updated for new directory structure**

### Julia Tests
- Smart path detection automatically finds `Inctools` package
- Work from both `julia/` and `tests/` directories
- Clear error messages if run from wrong location

### R Tests
- Source from `InctoolsJulia/R/` (updated from `R/`)
- Must run from `julia/` directory
- Work with reorganized file structure

### Documentation
- Clear usage instructions in `TEST_README.md`
- All changes documented in this file

**Status:** Ready to use from any supported location
