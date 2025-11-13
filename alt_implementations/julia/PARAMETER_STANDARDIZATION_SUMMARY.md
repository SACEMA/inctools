# Parameter Standardization: `cov` → `covar`

**Date:** 2025-11-13
**Status:** ✅ Complete and tested

---

## Summary

Successfully standardized the parameter naming across all `inccounts` functions from inconsistent `cov`/`covar` to consistent `covar`, while maintaining full backwards compatibility.

---

## Changes Made

### 1. Single Survey `inccounts()` (Lines 1155-1186)

**Before:**
```julia
function inccounts(n::Int64, npos::Int64, ntestR::Int64, nR::Int64,
    mdri::Float64, frr::Float64;
    # ... other parameters ...
    cov::Float64 = 0.0,  # ← Inconsistent naming
    # ... more parameters ...)

    # BUG: Hardcoded 0.0, ignored user input!
    return incprops(prev, prevR, mdri, frr, σ_prev = σ_prev,
                    σ_prevR = σ_prevR, σ_mdri = σ_mdri, σ_frr = σ_frr,
                    cov = 0.0,  # ← BUG!
                    T = T, timeconversion = timeconversion,
                    bs = bs, α = α, per = per)
end
```

**After:**
```julia
function inccounts(n::Int64, npos::Int64, ntestR::Int64, nR::Int64,
    mdri::Float64, frr::Float64;
    # ... other parameters ...
    covar::Union{Float64, Nothing} = nothing,  # ← New standard parameter
    cov::Union{Float64, Nothing} = nothing,    # ← Deprecated with warning
    # ... more parameters ...)

    # Deprecation handling logic
    if cov !== nothing && covar !== nothing
        if cov != covar
            error("Both 'cov' and 'covar' were provided with different values...")
        end
        @warn "The 'cov' argument is deprecated..."
        covar_value = covar
    elseif cov !== nothing
        @warn "The 'cov' argument is deprecated..."
        covar_value = cov
    elseif covar !== nothing
        covar_value = covar
    else
        covar_value = 0.0  # Default
    end

    # BUG FIXED: Now uses actual value
    return incprops(prev, σ_prev, prevR, σ_prevR, mdri, σ_mdri, frr, σ_frr,
                    covar = covar_value,  # ← Fixed!
                    T = T, timeconversion = timeconversion,
                    bs = bs, α = α, per = per)
end
```

### 2. Multiple Surveys `inccounts()` (Lines 1200-1232)

**Before:**
```julia
function inccounts(n::AbstractVector{Int64}, npos::AbstractVector{Int64},
    ntestR::AbstractVector{Int64}, nR::AbstractVector{Int64},
    mdri::Float64, frr::Float64;
    # ... other parameters ...
    # BUG: Referenced undefined `prev` variable!
    cov::Array{Float64,2} = Matrix{Float64}(I, size(prev)[1], size(prev)[1]),
    # ... more parameters ...)

    # BUG: Hardcoded 0.0!
    return incprops(prev, prevR, mdri, frr, σ_prev = σ_prev,
                    σ_prevR = σ_prevR, σ_mdri = σ_mdri, σ_frr = σ_frr,
                    cov = 0.0,  # ← BUG!
                    ...)
end
```

**After:**
```julia
function inccounts(n::AbstractVector{Int64}, npos::AbstractVector{Int64},
    ntestR::AbstractVector{Int64}, nR::AbstractVector{Int64},
    mdri::Float64, frr::Float64;
    # ... other parameters ...
    covar::Union{Array{Float64,2}, Nothing} = nothing,  # ← Fixed default
    cov::Union{Array{Float64,2}, Nothing} = nothing,    # ← Deprecated
    # ... more parameters ...)

    # Compute prev first
    prev, σ_prev  = prevalence(npos, n, de_npos)
    prevR, σ_prevR = prevalence(nR, ntestR, de_nR)

    # Deprecation handling (same logic as single survey)
    # ... (omitted for brevity) ...

    if covar === nothing && cov === nothing
        # BUG FIXED: Now computes after prev is defined
        covar_value = Matrix{Float64}(I, size(prev)[1], size(prev)[1])
    end

    # BUG FIXED: Uses actual value
    return incprops(prev, σ_prev, prevR, σ_prevR, mdri, σ_mdri, frr, σ_frr,
                    covar = covar_value,  # ← Fixed!
                    ...)
end
```

### 3. Additional Bug Fix: Corrected `incprops()` Call Signature

**Issue:** Both `inccounts` functions were calling `incprops` with incorrect parameter order.

**Before:**
```julia
incprops(prev, prevR, mdri, frr,   # ← Wrong order
         σ_prev = σ_prev,          # ← Wrong: keyword args
         σ_prevR = σ_prevR,
         σ_mdri = σ_mdri,
         σ_frr = σ_frr,
         ...)
```

**After:**
```julia
incprops(prev, σ_prev, prevR, σ_prevR, mdri, σ_mdri, frr, σ_frr,  # ← Correct!
         covar = covar_value,     # ← Keyword args start here
         ...)
```

**Reason:** The `incprops()` signature requires interleaved mean/σ positional arguments:
```julia
function incprops(prev::Float64,
    σ_prev::Float64,      # ← Positional
    prevR::Float64,
    σ_prevR::Float64,     # ← Positional
    mdri::Float64,
    σ_mdri::Float64,      # ← Positional
    frr::Float64,
    σ_frr::Float64;       # ← Semicolon: keyword args start below
    covar::Float64 = 0.0,
    ...)
```

---

## Bugs Fixed

### Bug 1: `cov` parameter ignored in single survey
- **Location:** Line 1168 (before fix)
- **Issue:** User could pass `cov=0.5` but it would be ignored; hardcoded `0.0` was used
- **Impact:** Correlation between prev and prevR was never applied in `inccounts`
- **Fix:** Now passes `covar_value` to `incprops`

### Bug 2: `cov` parameter ignored in multiple surveys
- **Location:** Line 1213 (before fix)
- **Issue:** Same as Bug 1 for multiple survey case
- **Fix:** Now passes `covar_value` to `incprops`

### Bug 3: Default value referenced undefined variable
- **Location:** Line 1200 (before fix)
- **Issue:** `cov::Array{Float64,2} = Matrix{Float64}(I, size(prev)[1], size(prev)[1])`
  referenced `prev` before it was computed
- **Impact:** Would cause an error if Julia ever evaluated the default expression
- **Fix:** Changed default to `nothing`, compute identity matrix after `prev` is available

### Bug 4: Incorrect `incprops()` call signature
- **Location:** Lines 1166-1169, 1211-1213 (before fix)
- **Issue:** Called `incprops` with means as positional args, σ's as keyword args
- **Impact:** Would fail with MethodError (somehow wasn't caught in previous testing?)
- **Fix:** Corrected to interleaved positional arguments

---

## Backwards Compatibility

✅ **Fully backwards compatible**

- Users can continue using `cov` parameter (with deprecation warning)
- New code should use `covar` parameter
- Both parameters cannot be provided with different values (error)
- Default behavior unchanged (0.0 for scalar, identity matrix for array)

### Deprecation Warning Example
```
┌ Warning: The 'cov' argument is deprecated and will be removed in the
│ next version. Please use 'covar' instead.
└ @ Inctools ~/dev/.../Inctools/src/Inctools.jl:1171
```

---

## Testing

Created comprehensive test suite: `test_covar_deprecation.jl`

### All 6 Tests Passed:

1. ✅ **New `covar` parameter** - Works correctly, no warning
2. ✅ **Deprecated `cov` parameter** - Works with deprecation warning
3. ✅ **Neither parameter** - Defaults to 0.0 correctly
4. ✅ **Both with same value** - Accepts with single warning
5. ✅ **Both with different values** - Correctly errors
6. ✅ **Bug fix verification** - Confidence intervals differ with different `covar` values

**Key result from Test 6:**
```
With covar=0.00015: I=0.0752, CI=[0.0343, 0.126]
With covar=0.0:     I=0.0752, CI=[0.0385, 0.123]
✓ Confidence intervals differ - covar is being used correctly!
```

This proves the bug is fixed - previously both would give identical CIs.

---

## Parameter Naming Now Consistent Across Package

| Function     | Parameter Name | Type              | Status       |
|--------------|----------------|-------------------|--------------|
| `incprops`   | `covar`        | Float64           | ✅ Standard  |
| `incdif`     | `covar`        | AbstractVector    | ✅ Standard  |
| `inccounts`  | `covar`        | Float64 or Array  | ✅ Standard  |
| `inccounts`  | `cov`          | Float64 or Array  | ⚠️ Deprecated |

---

## Migration Guide for Users

### Old Code (will still work with warning):
```julia
result = inccounts(1000, 200, 180, 20, 130.0, 0.01,
                   cov=0.0002, bs=1000)
```

### New Code (recommended):
```julia
result = inccounts(1000, 200, 180, 20, 130.0, 0.01,
                   covar=0.0002, bs=1000)
```

---

## Files Modified

1. **`Inctools/src/Inctools.jl`**
   - Lines 1155-1186: Single survey `inccounts()` function
   - Lines 1200-1232: Multiple survey `inccounts()` function

2. **Tests Created:**
   - `test_covar_deprecation.jl` - Comprehensive deprecation testing

---

## Conclusion

✅ **All objectives achieved:**
- Parameter naming standardized to `covar` across entire package
- Full backwards compatibility maintained with deprecation warnings
- **4 significant bugs fixed** (ignored parameter values + incorrect function calls)
- Comprehensive test coverage confirms correct behavior
- All existing tests continue to pass

**Status:** Ready for production use
