# R Package Renaming: Inctools → InctoolsJulia

**Date:** 2025-11-13
**Reason:** Avoid conflict with existing 'inctools' package on CRAN
**Status:** ✅ Complete

---

## Issue

The original R package wrapper was named **Inctools**, which could conflict with the existing **inctools** package on CRAN:
- While technically different (capital I vs lowercase i), this would be confusing for users
- Could cause conflicts or unexpected behavior
- Users might accidentally load the wrong package

## Solution

Renamed the R package to **InctoolsJulia** throughout:
- Clearly distinguishes it from CRAN's inctools package
- Makes it obvious this is the Julia backend version
- Allows users to have both packages installed without conflicts
- Follows R package naming conventions

---

## Files Updated

### 1. DESCRIPTION
- **Before:** `Package: Inctools`
- **After:** `Package: InctoolsJulia`
- Added note in description: "distinct from the 'inctools' R package on CRAN"

### 2. R/zzz.R
- **Before:** `package = "Inctools"`
- **After:** `package = "InctoolsJulia"`
- Updated path detection for installed package

### 3. README_R.md
- **Before:** `library(Inctools)`
- **After:** `library(InctoolsJulia)`
- Updated title to "InctoolsJulia"
- Added prominent note about naming distinction
- Updated all usage examples

### 4. R_API_SUMMARY.md
- Updated package name references
- Clarified comparison table
- Updated conclusion section
- Changed all `library(Inctools)` to `library(InctoolsJulia)`

### 5. install_R_package.R
- Updated title comment
- Updated startup message

### 6. R/inctools.R
- No changes needed (function names remain the same)

### 7. NAMESPACE
- No changes needed (exported functions remain the same)

### 8. test_R_api.R
- No changes needed (uses source(), not library())

---

## User Impact

### Before (Confusing)
```r
# Which inctools package is this?
library(Inctools)  # Julia version
library(inctools)  # CRAN version - CONFLICT!
```

### After (Clear)
```r
# Clear distinction
library(InctoolsJulia)  # Julia backend version
library(inctools)       # CRAN version - no conflict
```

### Users can now:
1. Install both packages without conflicts
2. Easily distinguish which version they're using
3. Switch between implementations as needed

---

## Function Names

**No changes to function names** - all remain the same:
- `inctools_setup()`
- `prevalence()`
- `incprops()`
- `inccounts()`
- `incdif()`

Only the **package name** changed (what you put in `library()`).

---

## Installation Methods

### Method 1: Source files
```r
source("R/zzz.R")
source("R/inctools.R")
# No package name needed
```

### Method 2: Install as package
```r
install.packages(".", repos = NULL, type = "source")
library(InctoolsJulia)  # Changed from library(Inctools)
```

### Method 3: Installation script
```bash
Rscript install_R_package.R
# Will show: "InctoolsJulia - R Interface to Inctools.jl"
```

---

## Documentation Updates

All documentation now clearly states:
- R package name: **InctoolsJulia**
- Julia package name: **Inctools.jl** (unchanged)
- CRAN package name: **inctools** (external, unchanged)

### Example from README_R.md:
```
# InctoolsJulia - R Interface to Inctools.jl

**Note:** This package is named **InctoolsJulia** to distinguish it from the
existing **inctools** package on CRAN. Both packages provide similar functionality,
but this version uses a high-performance Julia backend.
```

---

## Comparison Table

| Package | Source | Backend | Installation | Name in R |
|---------|--------|---------|--------------|-----------|
| **inctools** | CRAN | Pure R | `install.packages("inctools")` | `library(inctools)` |
| **InctoolsJulia** | This repo | Julia | Manual install | `library(InctoolsJulia)` |
| **Inctools.jl** | Julia | Julia | Julia Pkg manager | `using Inctools` |

---

## Testing

All functionality remains identical:
- ✅ Function names unchanged
- ✅ Function signatures unchanged
- ✅ API behavior unchanged
- ✅ Only the package loading name changed

Test with:
```r
Rscript test_R_api.R
# All tests still pass
```

---

## Benefits of New Naming

1. **No conflicts** - Can coexist with CRAN's inctools
2. **Clear identity** - Obvious it's the Julia version
3. **Professional** - Follows R package naming best practices
4. **Flexible** - Users can choose which to use
5. **Safe** - No risk of loading wrong package

---

## Checklist

- [x] Updated DESCRIPTION (package name)
- [x] Updated R/zzz.R (package path)
- [x] Updated README_R.md (all references)
- [x] Updated R_API_SUMMARY.md (all references)
- [x] Updated install_R_package.R (messages)
- [x] Verified test_R_api.R (no changes needed)
- [x] Verified R/inctools.R (no changes needed)
- [x] Verified NAMESPACE (no changes needed)
- [x] Created this summary document

---

## Conclusion

✅ **Package successfully renamed from Inctools to InctoolsJulia**

The renaming:
- Eliminates potential conflicts with CRAN's inctools package
- Makes the package identity clear and unambiguous
- Maintains all functionality (only loading name changed)
- Allows users to install both packages if desired

**No further action needed** - all files updated and consistent.
