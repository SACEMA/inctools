# Directory Reorganization Summary

**Date:** 2025-11-13
**Status:** ✅ Complete

---

## Motivation

The original directory structure mixed R and Julia files at the top level, making it:
- Confusing which files belonged to which package
- Difficult to distribute packages separately
- Non-standard for R package conventions
- Messy with test files scattered around

---

## Changes Made

### Before (Messy)

```
julia/
├── Inctools/              # Julia package (subdirectory)
│   ├── Project.toml
│   └── src/
├── R/                     # R files at top level
├── DESCRIPTION            # R files at top level
├── NAMESPACE
├── README_R.md
├── install_R_package.R
├── test_comprehensive.jl  # Tests scattered
├── test_R_api.R
├── test_rtmvnorm.jl
└── ...
```

### After (Clean)

```
julia/
├── Inctools/              # Julia package
│   ├── Project.toml
│   ├── Manifest.toml
│   └── src/
│       └── Inctools.jl
│
├── InctoolsJulia/         # R package
│   ├── R/
│   │   ├── zzz.R
│   │   └── inctools.R
│   ├── inst/
│   │   └── Inctools -> ../../Inctools  # Symlink
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   ├── README.md
│   ├── install_R_package.R
│   ├── R_API_SUMMARY.md
│   ├── R_API_PARAMETER_FIX.md
│   └── PACKAGE_RENAME_SUMMARY.md
│
├── tests/                 # All tests together
│   ├── test_comprehensive.jl
│   ├── test_R_api.R
│   ├── test_rtmvnorm.jl
│   ├── test_gibbs_vs_rejection.jl
│   ├── test_covar_deprecation.jl
│   └── test_R_fix.R
│
├── README.md              # Main README
├── REVIEW_SUMMARY.md
├── PARAMETER_STANDARDIZATION_SUMMARY.md
└── test.jl
```

---

## File Movements

### Moved to `InctoolsJulia/`
- ✅ `R/` → `InctoolsJulia/R/`
- ✅ `DESCRIPTION` → `InctoolsJulia/DESCRIPTION`
- ✅ `NAMESPACE` → `InctoolsJulia/NAMESPACE`
- ✅ `README_R.md` → `InctoolsJulia/README.md` (renamed)
- ✅ `install_R_package.R` → `InctoolsJulia/install_R_package.R`
- ✅ `R_API_*.md` → `InctoolsJulia/R_API_*.md`
- ✅ `PACKAGE_RENAME_SUMMARY.md` → `InctoolsJulia/PACKAGE_RENAME_SUMMARY.md`

### Moved to `tests/`
- ✅ `test_comprehensive.jl` → `tests/test_comprehensive.jl`
- ✅ `test_R_api.R` → `tests/test_R_api.R`
- ✅ `test_rtmvnorm.jl` → `tests/test_rtmvnorm.jl`
- ✅ `test_gibbs_vs_rejection.jl` → `tests/test_gibbs_vs_rejection.jl`
- ✅ `test_covar_deprecation.jl` → `tests/test_covar_deprecation.jl`
- ✅ `test_R_fix.R` → `tests/test_R_fix.R`

### Created New
- ✅ `InctoolsJulia/inst/` directory
- ✅ `InctoolsJulia/inst/Inctools` symlink → `../../Inctools`
- ✅ `README.md` (main README at top level)
- ✅ `tests/` directory

### Remained at Top Level
- ✅ `Inctools/` (Julia package)
- ✅ `test.jl` (simple Julia test)
- ✅ `REVIEW_SUMMARY.md`
- ✅ `PARAMETER_STANDARDIZATION_SUMMARY.md`

---

## Code Updates

### 1. `InctoolsJulia/R/zzz.R`

**Updated path detection logic:**

```r
# When installed: inst/Inctools becomes <package-root>/Inctools
pkg_path <- system.file("Inctools", package = "InctoolsJulia", mustWork = FALSE)

if (pkg_path == "") {
  # If not installed, try common development locations
  possible_paths <- c(
    file.path(getwd(), "..", "Inctools"),        # From InctoolsJulia/
    file.path(getwd(), "Inctools"),              # From julia/
    file.path(dirname(getwd()), "Inctools")      # From subdirectory
  )
  # ... search logic
}
```

Now works from multiple locations:
- When sourcing from `InctoolsJulia/` directory
- When running from `julia/` directory
- When installed as R package

### 2. `InctoolsJulia/install_R_package.R`

**Updated to find Inctools in multiple locations:**

```r
possible_paths <- c(
  file.path(getwd(), "..", "Inctools"),   # From InctoolsJulia/
  file.path(getwd(), "Inctools"),         # From julia/
  file.path(dirname(getwd()), "Inctools") # From subdirectory
)
```

**Updated installation instructions:**
- Run from `julia/` or `julia/InctoolsJulia/`
- Install command: `install.packages('InctoolsJulia', ...)`
- Test command: `Rscript tests/test_R_api.R`

---

## R Package Installation (inst/ Directory)

### How It Works

R packages use the `inst/` directory for files that should be installed with the package. During installation:

1. **Development:** `InctoolsJulia/inst/Inctools` is a symlink to `../../Inctools`
2. **Installation:** R copies `inst/` contents to the package root
3. **Installed:** Becomes `<R-library>/InctoolsJulia/Inctools/`

### Why Symlink?

- ✅ Avoids duplicating the entire Inctools.jl package
- ✅ Keeps Julia package as single source of truth
- ✅ Changes to Inctools.jl are automatically reflected
- ✅ Saves disk space

### Installation Behavior

```bash
# Before installation (development)
InctoolsJulia/
├── inst/
│   └── Inctools -> ../../Inctools  # Symlink follows to real files

# After installation
/Library/R/library/InctoolsJulia/
├── Inctools/                        # Real files copied
│   ├── Project.toml
│   └── src/
└── R/
```

---

## Usage Changes

### Julia (No Change)

```bash
# From julia/ directory
julia
julia> using Pkg
julia> Pkg.activate("./Inctools")
julia> using Inctools
```

### R - Option 1: Source Files

**Before:**
```r
# From julia/ directory
source("R/zzz.R")
source("R/inctools.R")
```

**After:**
```r
# From julia/InctoolsJulia/ directory
source("R/zzz.R")
source("R/inctools.R")
```

### R - Option 2: Install Package

**Before:**
```r
# From julia/ directory
install.packages(".", repos = NULL, type = "source")
library(Inctools)  # Wrong name
```

**After:**
```r
# From julia/ directory
install.packages("InctoolsJulia", repos = NULL, type = "source")
library(InctoolsJulia)  # Correct name
```

### Running Tests

**Before:**
```bash
julia test_comprehensive.jl
Rscript test_R_api.R
```

**After:**
```bash
julia tests/test_comprehensive.jl
Rscript tests/test_R_api.R
```

---

## Benefits

### 1. Clear Separation
- ✅ Julia package in `Inctools/`
- ✅ R package in `InctoolsJulia/`
- ✅ Tests in `tests/`
- ✅ Docs at appropriate levels

### 2. R Package Standards
- ✅ Follows R package conventions
- ✅ Uses `inst/` for bundled files
- ✅ Self-contained package directory
- ✅ Can be installed directly

### 3. Easy Distribution
- ✅ Each package can be zipped/distributed separately
- ✅ Clear which files belong to which package
- ✅ Simpler CI/CD setup

### 4. Better Organization
- ✅ All tests in one place
- ✅ Documentation clearly labeled
- ✅ No confusion about file ownership

---

## Directory Listing

```
julia/
├── Inctools/                          # Julia package
│   ├── Project.toml                   # v0.2.0
│   ├── Manifest.toml
│   └── src/
│       └── Inctools.jl                # 1234 lines
│
├── InctoolsJulia/                     # R package
│   ├── R/
│   │   ├── zzz.R                      # Package initialization
│   │   └── inctools.R                 # Main functions (219 lines)
│   ├── inst/
│   │   └── Inctools -> ../../Inctools # Symlink to Julia package
│   ├── DESCRIPTION                    # v0.2.0
│   ├── NAMESPACE
│   ├── README.md                      # R package docs
│   ├── install_R_package.R
│   ├── R_API_SUMMARY.md
│   ├── R_API_PARAMETER_FIX.md
│   └── PACKAGE_RENAME_SUMMARY.md
│
├── tests/                             # All tests
│   ├── test_comprehensive.jl          # 9 Julia tests
│   ├── test_rtmvnorm.jl
│   ├── test_gibbs_vs_rejection.jl
│   ├── test_covar_deprecation.jl
│   ├── test_R_api.R                   # 6 R tests
│   └── test_R_fix.R
│
├── README.md                          # Main README
├── REVIEW_SUMMARY.md                  # Code review
├── PARAMETER_STANDARDIZATION_SUMMARY.md
├── DIRECTORY_REORGANIZATION.md        # This file
└── test.jl                            # Simple test
```

---

## Verification

### Check Structure
```bash
ls -la
# Should see: Inctools/ InctoolsJulia/ tests/ README.md
```

### Verify Symlink
```bash
ls -la InctoolsJulia/inst/
# Should see: Inctools -> ../../Inctools
```

### Test Julia Package
```bash
julia tests/test_comprehensive.jl
# Should pass all 9 tests
```

### Test R Package
```bash
Rscript tests/test_R_api.R
# Should pass all 6 tests
```

---

## Migration Checklist

- [x] Create `InctoolsJulia/` directory
- [x] Create `tests/` directory
- [x] Move R files to `InctoolsJulia/`
- [x] Move test files to `tests/`
- [x] Create `InctoolsJulia/inst/` directory
- [x] Create symlink `inst/Inctools -> ../../Inctools`
- [x] Update `R/zzz.R` path detection
- [x] Update `install_R_package.R` paths
- [x] Rename `README_R.md` → `README.md`
- [x] Create main `README.md`
- [x] Verify all tests still pass
- [x] Document reorganization

---

## Conclusion

✅ **Directory reorganization complete!**

The new structure:
- Clearly separates Julia and R packages
- Follows R package conventions (inst/ directory)
- Organizes tests in dedicated directory
- Makes distribution easier
- Maintains all functionality
- All tests pass

**No code changes needed** - only directory structure and path references updated.

**Status:** Ready for use in new structure
