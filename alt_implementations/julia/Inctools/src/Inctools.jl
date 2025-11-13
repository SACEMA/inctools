# Incidence Estimation Tools (Julia implementation). Copyright (C) 2019,
# Eduard Grebe and individual contributors.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

module Inctools

export prevalence, rtmvnorm, incprops, inccounts, incdif

import LinearAlgebra: I, diag
import Distributions
import Statistics
import DataFrames

"""
    prevalence(pos, n, de=1; ci=false, α=0.05)

Estimate prevalence from the number of positive results in a sample.

This function calculates point estimates and standard errors for prevalence using the
normal approximation. Optionally provides confidence intervals using exact F-distribution
methods (Clopper-Pearson) for boundary cases.

# Arguments
- `pos::Int`: Number of positive test results
- `n::Int`: Total sample size
- `de::Float64=1`: Design effect for complex survey designs (default: 1 for simple random sampling)
- `ci::Bool=false`: Whether to calculate confidence intervals
- `α::Float64=0.05`: Significance level for confidence intervals (default: 0.05 for 95% CI)

# Returns
A named tuple containing:
- `p::Float64`: Estimated prevalence (proportion)
- `σ::Float64`: Standard error of the prevalence estimate
- `ci::Vector{Float64}`: Confidence interval [lower, upper] (only if `ci=true`)

# Notes
- Issues a warning if `pos < 5` or `n - pos < 5` as the normal approximation may not be valid
- Uses exact F-distribution method (Clopper-Pearson) for confidence intervals when `pos = 0` or `pos = n`
- If `n = 0`, automatically sets `n = 1` and issues a warning

# Examples
```julia
# Simple prevalence estimate
prevalence(150, 1000)
# Output: (p = 0.15, σ = 0.0112915...)

# With confidence intervals
prevalence(150, 1000, ci=true)
# Output: (p = 0.15, σ = 0.0112915..., ci = [0.128..., 0.174...])

# Accounting for design effect
prevalence(150, 1000, 2.0, ci=true)
# Output: (p = 0.15, σ = 0.0159..., ci = [0.128..., 0.174...])
```
"""
function prevalence(pos, n, de = 1; ci = false, α = 0.05) #, f = 1
    if n == 0
        @warn "n = 0: prevalence undefined. n set to 1."
        n = 1
    end
    if pos < 5
        @warn "Too few successes for the normal approximation to be valid"
    end
    if n - pos < 5
        @warn "Sample size too small for the normal approximation to be valid"
    end
    p = pos/n
    σ = sqrt( (p * (1 - p)) / n ) * de #* sqrt(1 - f)
    if !ci
        return (p = p, σ = σ)
    elseif ci
        Fd =
        if pos == 0
            lb = 0
            ub = 1-(α/2)^(1/n)
        elseif pos == n
            lb = (α/2)^(1/n)
            ub = 1
        else
            lb = ( 1+(n-pos+1)/(pos * Distributions.quantile(Distributions.FDist(2*pos, 2*(n-pos+1)), α/2)) ) ^ (-1)
            ub = ( 1+(n-pos)/((pos + 1) * Distributions.quantile(Distributions.FDist(2*(pos + 1), 2*(n-pos)), 1-α/2)) ) ^ (-1)
        end

        # Normal approximation
        # pprime = (pos + 2) / (n + 4)
        # sigmaprime = sqrt( (pprime * (1 - pprime)) / (n + 4) )
        #lb = max(pprime - 1.96 * sigmaprime,0)
        # ub = min(pprime + 1.96 * sigmaprime,1)
        return (p = p, σ = σ, ci = [lb, ub])
    end
end

function kassanjee(prev::Float64, prevR::Float64, mdri::Float64, frr::Float64, T::Float64)
    return (prev * (prevR - frr)) / ((1 - prev) * (mdri - frr * T))
end

function σ_dm(prev, prevR, mdri, frr, T, σ_prev, σ_prevR, σ_mdri, σ_frr)
    fot_prev = (prevR - frr) / (((1 - prev)^2) * (mdri - frr * T))
    fot_prevR = prev / ((1 - prev) * (mdri - frr * T))
    fot_mdri = (frr * prev - prevR * prev) / ((1 - prev) * ((mdri - frr * T)^2))
    fot_frr = (prev * (T * prevR - mdri)) / ((1 - prev) * ((mdri - frr * T)^2))
    σ = sqrt(fot_prev^2 * σ_prev^2 + fot_prevR^2 * σ_prevR^2 + fot_mdri^2 * σ_mdri^2 + fot_frr^2 * σ_frr^2)
    σ_infSS = sqrt(fot_mdri^2 * σ_mdri^2 + fot_frr^2 * σ_frr^2)
    return σ, σ_infSS
end

function σ_Δ_dm(prev, prevR, mdri, frr, T, σ_prev, σ_prevR, σ_mdri, σ_frr)
    fot_prev1 = (prevR[1] - frr) / (((1 - prev[1])^2) * (mdri - frr * T))
    fot_prevR1 = prev[1] / ((1 - prev[1]) * (mdri - frr * T))
    fot_mdri1 = (frr * prev[1] - prevR[1] * prev[1]) / ((1 - prev[1]) * ((mdri - frr * T)^2))
    fot_frr1 = (prev[1] * (T * prevR[1] - mdri)) / ((1 - prev[1]) * ((mdri - frr * T)^2))

    fot_prev2 = (prevR[2] - frr) / (((1 - prev[2])^2) * (mdri - frr * T))
    fot_prevR2 = prev[2] / ((1 - prev[2]) * (mdri - frr * T))
    fot_mdri2 = (frr * prev[2] - prevR[2] * prev[2]) / ((1 - prev[2]) * ((mdri - frr * T)^2))
    fot_frr2 = (prev[2] * (T * prevR[2] - mdri)) / ((1 - prev[2]) * ((mdri - frr * T)^2))

    variance = ((fot_prev1^2) * σ_prev[1]^2) + ((fot_prev2^2) * σ_prev[2]^2) + ((fot_prevR1^2) * σ_prevR[1]^2) + ((fot_prevR2^2) * σ_prevR[2]^2) + ((fot_mdri1 - fot_mdri2)^2 * σ_mdri^2) + ((fot_frr1 - fot_frr2)^2 * σ_frr^2)
    σ = sqrt(variance)
    return σ
end

# This follows the logic of Stefan Wilhelm
# https://github.com/cran/tmvtnorm/blob/master/R/rtmvnorm.R
function rtnorm_gibbs(n::Int64, μ::Float64, σ::Float64, lower::Float64, upper::Float64)
    d = Distributions.Normal(μ,σ)
    F = rand(Distributions.Uniform(), n)
    lcdf = Distributions.cdf(d, lower)
    ucdf = Distributions.cdf(d, upper)
    tp = ucdf - lcdf
    Q = F .* tp  .+ lcdf
    r = Distributions.quantile.(Distributions.Normal(0, 1), Q) .* σ .+ μ
    return r
end

"""
    rtmvnorm_rejection(n, µ, Σ, lower, upper)

Generate random samples from a truncated multivariate normal distribution using rejection sampling.

⚠️ **Warning**: This function is hardcoded for 4-dimensional distributions. It uses rejection
sampling which becomes very inefficient when the truncation region has low probability mass.
For diagonal covariance matrices, consider using `rtmvnorm_gibbs` instead.

# Arguments
- `n::Int64`: Number of samples to generate
- `µ::AbstractVector{Float64}`: Mean vector (length 4)
- `Σ::Array{Float64,2}`: Covariance matrix (4×4)
- `lower::AbstractVector{Float64}`: Lower truncation bounds (length 4)
- `upper::AbstractVector{Float64}`: Upper truncation bounds (length 4)

# Returns
- `Matrix{Float64}`: n×4 matrix of samples, where each row is one sample

# Algorithm
Uses rejection sampling:
1. Estimates acceptance rate with 1000 trial samples
2. Draws samples from untruncated MVN
3. Rejects samples outside truncation bounds
4. Repeats until n valid samples obtained

# Performance Notes
- Becomes very inefficient when truncation region has low probability
- Acceptance rate estimation adds overhead but reduces wasted sampling
- Consider `rtmvnorm_gibbs` for diagonal covariance matrices

# Examples
```julia
using LinearAlgebra

# Sample from 4D normal with correlation
µ = [0.5, 0.5, 0.5, 0.5]
Σ = [0.1 0.02 0 0; 0.02 0.1 0 0; 0 0 0.1 0; 0 0 0 0.1]
samples = rtmvnorm_rejection(1000, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])
```

# See Also
- `rtmvnorm`: Smart wrapper that chooses between Gibbs and rejection sampling
- `rtmvnorm_gibbs`: Efficient Gibbs sampler for diagonal covariance
- `rtnorm_gibbs`: Univariate Gibbs sampler
"""
# Rejection sampling implementation - retained for correlated variables
function rtmvnorm_rejection(n::Int64,
    µ::AbstractVector{Float64},
    Σ::Array{Float64,2},
    lower::AbstractVector{Float64},
    upper::AbstractVector{Float64})

    d = Distributions.MvNormal(µ, Σ)

    # find acceptance rate
    r = transpose(rand(d, 1000))
    racc = r[(r[:,1] .>= lower[1]) .&
             (r[:,2] .>= lower[2]) .&
             (r[:,3] .>= lower[3]) .&
             (r[:,4] .>= lower[4]) .&
             (r[:,1] .<= upper[1]) .&
             (r[:,2] .<= upper[2]) .&
             (r[:,3] .<= upper[3]) .&
             (r[:,4] .<= upper[4]),:]
    rr = 1 - size(racc,1)/1000

    r = transpose(rand(d, Int(round(n + n * rr))))
    r = r[(r[:,1] .>= lower[1]) .&
             (r[:,2] .>= lower[2]) .&
             (r[:,3] .>= lower[3]) .&
             (r[:,4] .>= lower[4]) .&
             (r[:,1] .<= upper[1]) .&
             (r[:,2] .<= upper[2]) .&
             (r[:,3] .<= upper[3]) .&
             (r[:,4] .<= upper[4]),:]
    accepted = size(r,1)

    while accepted < n
        radd = transpose(rand(d, Int(round(n * rr + 1))))
        radd = radd[(radd[:,1] .>= lower[1]) .&
                    (radd[:,2] .>= lower[2]) .&
                    (radd[:,3] .>= lower[3]) .&
                    (radd[:,4] .>= lower[4]) .&
                    (radd[:,1] .<= upper[1]) .&
                    (radd[:,2] .<= upper[2]) .&
                    (radd[:,3] .<= upper[3]) .&
                    (radd[:,4] .<= upper[4]),:]
        r = [r ; radd]
        accepted = size(r,1)
    end
    r = r[1:n,:]
    return r
end

"""
    rtmvnorm_gibbs(n, µ, Σ, lower, upper)

Generate random samples from a truncated multivariate normal using independent Gibbs sampling.

This function uses Gibbs sampling when the covariance matrix is diagonal (no correlation between
variables). It samples each dimension independently using the efficient univariate `rtnorm_gibbs`
method. Much faster than rejection sampling for diagonal covariance matrices.

# Arguments
- `n::Int64`: Number of samples to generate
- `µ::AbstractVector{Float64}`: Mean vector (any length)
- `Σ::Array{Float64,2}`: Diagonal covariance matrix (d×d)
- `lower::AbstractVector{Float64}`: Lower truncation bounds (length d)
- `upper::AbstractVector{Float64}`: Upper truncation bounds (length d)

# Returns
- `Matrix{Float64}`: n×d matrix of samples, where each row is one sample

# Requirements
- The covariance matrix must be diagonal (off-diagonal elements = 0)
- Dimensions of µ, Σ, lower, and upper must match

# Algorithm
For each dimension i:
1. Extract µᵢ and σᵢ = √Σᵢᵢ
2. Use univariate Gibbs sampler: `rtnorm_gibbs(n, µᵢ, σᵢ, lowerᵢ, upperᵢ)`
3. Combine samples into matrix

# Performance
- O(n × d) complexity
- Much faster than rejection sampling for diagonal covariance
- No wasted samples (100% acceptance rate)

# Examples
```julia
using LinearAlgebra

# 4D independent truncated normals
µ = [0.5, 0.3, 0.7, 0.4]
Σ = Matrix(Diagonal([0.1, 0.08, 0.12, 0.09]))
samples = rtmvnorm_gibbs(1000, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])

# 2D case
µ = [0.2, 0.8]
Σ = Matrix(Diagonal([0.05, 0.05]))
samples = rtmvnorm_gibbs(5000, µ, Σ, [0.0, 0.0], [1.0, 1.0])
```

# See Also
- `rtmvnorm`: Smart wrapper that automatically chooses this method for diagonal covariance
- `rtmvnorm_rejection`: For correlated variables (non-diagonal covariance)
- `rtnorm_gibbs`: Univariate Gibbs sampler used internally
"""
function rtmvnorm_gibbs(n::Int64,
    µ::AbstractVector{Float64},
    Σ::Array{Float64,2},
    lower::AbstractVector{Float64},
    upper::AbstractVector{Float64})

    d = length(µ)

    # Verify dimensions match
    if size(Σ) != (d, d)
        throw(DimensionMismatch("Covariance matrix dimensions $(size(Σ)) do not match mean vector length $d"))
    end
    if length(lower) != d || length(upper) != d
        throw(DimensionMismatch("Truncation bounds must have same length as mean vector ($d)"))
    end

    # Extract standard deviations from diagonal
    σ = sqrt.(diag(Σ))

    # Sample each dimension independently using Gibbs sampler
    samples = zeros(n, d)
    for i in 1:d
        samples[:, i] = rtnorm_gibbs(n, µ[i], σ[i], lower[i], upper[i])
    end

    return samples
end

"""
    rtmvnorm(n, µ, Σ, lower, upper; method=:auto)

Generate random samples from a truncated multivariate normal distribution.

This is a smart wrapper that automatically chooses between Gibbs sampling (for diagonal covariance)
and rejection sampling (for correlated variables). Defaults to Gibbs sampling when possible, as it
is much more efficient.

# Arguments
- `n::Int64`: Number of samples to generate
- `µ::AbstractVector{Float64}`: Mean vector (length d)
- `Σ::Array{Float64,2}`: Covariance matrix (d×d)
- `lower::AbstractVector{Float64}`: Lower truncation bounds (length d)
- `upper::AbstractVector{Float64}`: Upper truncation bounds (length d)
- `method::Symbol=:auto`: Sampling method - `:auto`, `:gibbs`, or `:rejection`
  - `:auto` - Automatically choose Gibbs for diagonal Σ, rejection otherwise
  - `:gibbs` - Force Gibbs sampling (error if Σ not diagonal)
  - `:rejection` - Force rejection sampling (works but slow for 4D only)

# Returns
- `Matrix{Float64}`: n×d matrix of samples, where each row is one sample

# Method Selection (when method=:auto)
- **Gibbs sampling** used when covariance matrix is diagonal (max off-diagonal < 1e-10)
- **Rejection sampling** used when variables are correlated (non-diagonal covariance)

# Notes
- Gibbs sampling: Fast, works for any dimension with diagonal covariance
- Rejection sampling: Slow, currently hardcoded for 4 dimensions only
- Gibbs is typically 10-100x faster than rejection for diagonal covariance

# Examples
```julia
using LinearAlgebra

# Example 1: Diagonal covariance (automatically uses Gibbs)
µ = [0.5, 0.3, 0.7, 0.4]
Σ = Matrix(Diagonal([0.1, 0.08, 0.12, 0.09]))
samples = rtmvnorm(1000, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])

# Example 2: Correlated variables (automatically uses rejection)
Σ_corr = [0.1 0.02 0 0; 0.02 0.08 0 0; 0 0 0.12 0; 0 0 0 0.09]
samples = rtmvnorm(1000, µ, Σ_corr, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])

# Example 3: Force specific method
samples = rtmvnorm(1000, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0], method=:gibbs)
```

# See Also
- `rtmvnorm_gibbs`: Efficient Gibbs sampler for diagonal covariance
- `rtmvnorm_rejection`: Rejection sampler for correlated variables
- `rtnorm_gibbs`: Univariate Gibbs sampler
"""
function rtmvnorm(n::Int64,
    µ::AbstractVector{Float64},
    Σ::Array{Float64,2},
    lower::AbstractVector{Float64},
    upper::AbstractVector{Float64};
    method::Symbol = :auto)

    d = length(µ)

    # Check if covariance matrix is diagonal (no correlation)
    is_diagonal = true
    for i in 1:d
        for j in 1:d
            if i != j && abs(Σ[i,j]) > 1e-10
                is_diagonal = false
                break
            end
        end
        if !is_diagonal
            break
        end
    end

    # Choose method
    if method == :auto
        if is_diagonal
            # Use efficient Gibbs sampling
            return rtmvnorm_gibbs(n, µ, Σ, lower, upper)
        else
            # Use rejection sampling for correlated variables
            if d != 4
                @warn "Rejection sampling currently only supports 4 dimensions. Covariance is non-diagonal, falling back to rejection sampling which requires exactly 4 dimensions."
            end
            return rtmvnorm_rejection(n, µ, Σ, lower, upper)
        end
    elseif method == :gibbs
        if !is_diagonal
            throw(ArgumentError("Gibbs sampling requires diagonal covariance matrix. Use method=:rejection or method=:auto"))
        end
        return rtmvnorm_gibbs(n, µ, Σ, lower, upper)
    elseif method == :rejection
        if d != 4
            throw(DimensionMismatch("Rejection sampling currently only supports 4 dimensions, got $d"))
        end
        return rtmvnorm_rejection(n, µ, Σ, lower, upper)
    else
        throw(ArgumentError("Unknown method: $method. Use :auto, :gibbs, or :rejection"))
    end
end

"""
    incprops(prev, σ_prev, prevR, σ_prevR, mdri, σ_mdri, frr, σ_frr; kwargs...)

Estimate HIV incidence from survey proportions using the Kassanjee method.

This function implements cross-sectional incidence estimation using prevalence of infection,
prevalence of recent infection (recency), mean duration of recent infection (MDRI), and
false-recent rate (FRR). Supports both Delta method and bootstrap approaches for
uncertainty estimation.

# Arguments
## Required Parameters
- `prev::Float64`: Prevalence of HIV infection (proportion, 0-1)
- `σ_prev::Float64`: Standard error of prevalence estimate
- `prevR::Float64`: Prevalence of recency among HIV-positive individuals (proportion, 0-1)
- `σ_prevR::Float64`: Standard error of recency prevalence
- `mdri::Float64`: Mean duration of recent infection (in same units as `T`, typically days)
- `σ_mdri::Float64`: Standard error of MDRI estimate
- `frr::Float64`: False-recent rate (proportion, 0-1)
- `σ_frr::Float64`: Standard error of FRR estimate

## Keyword Arguments
- `covar::Float64=0.0`: Covariance between prev and prevR (typically 0 for independent estimates)
- `T::Float64=730.5`: Time cutoff for recency test (default 730.5 days = 2 years)
- `timeconversion::Float64=365.25`: Conversion factor to express incidence per unit time
  (default 365.25 converts from daily to annual incidence)
- `bs::Int64=0`: Number of bootstrap iterations (0 = use Delta method only)
- `gibbs::Bool=false`: Use Gibbs sampler for truncated normals (false = use Distributions.jl)
- `bs_numbers::Bool=false`: Bootstrap from counts rather than proportions
- `bs_numbers_n::Vector{Int64}=[0,0]`: Sample sizes [n_total, n_positive] for count-based bootstrap
- `α::Float64=0.05`: Significance level for confidence intervals (default 0.05 = 95% CI)
- `per::Int64=1`: Multiplier for incidence rate (e.g., per=100 for percentage)

# Returns
A named tuple containing:
- `I::Float64`: Point estimate of incidence
- `CI::Vector{Float64}`: Confidence interval [lower, upper]
- `σ::Float64`: Standard error of incidence estimate
- `RSE::Float64`: Relative standard error (σ/|I|)
- `cov_prev_I::Float64`: Covariance between prevalence and incidence (bootstrap only)
- `cor_prev_I::Float64`: Correlation between prevalence and incidence (bootstrap only)

# Methods
1. **Delta Method** (`bs=0`): Fast analytical variance using first-order Taylor expansion
2. **Bootstrap** (`bs>0`): Resampling-based uncertainty estimation with truncated normals
3. **Count Bootstrap** (`bs_numbers=true`): Bootstrap from binomial counts for finite sample sizes

# Mathematical Model
The Kassanjee estimator is:
```
I = P × (P_R - FRR) / [(1 - P) × (MDRI - FRR × T)]
```
where P = prevalence, P_R = recency prevalence among positives

# Examples
```julia
# Example 1: Simple incidence estimate using Delta method
result = incprops(
    0.20, 0.015,      # Prevalence = 20% ± 1.5%
    0.10, 0.02,       # Recency = 10% ± 2%
    130.0, 15.0,      # MDRI = 130 days ± 15
    0.01, 0.005       # FRR = 1% ± 0.5%
)
println("Incidence: \$(result.I) per year (95% CI: \$(result.CI))")

# Example 2: Bootstrap with 10,000 iterations
result_bs = incprops(
    0.20, 0.015, 0.10, 0.02, 130.0, 15.0, 0.01, 0.005,
    bs = 10000
)

# Example 3: Account for covariance between prevalence estimates
result_cov = incprops(
    0.20, 0.015, 0.10, 0.02, 130.0, 15.0, 0.01, 0.005,
    covar = 0.0002,  # Positive covariance
    bs = 10000
)

# Example 4: Express incidence as percentage per 100 person-years
result_pct = incprops(
    0.20, 0.015, 0.10, 0.02, 130.0, 15.0, 0.01, 0.005,
    per = 100
)
```

# References
Kassanjee R, et al. (2012). A new general biomarker-based incidence estimator.
Epidemiology, 23(5), 721-728.

# See Also
- `inccounts`: Estimate incidence from counts rather than proportions
- `incdif`: Estimate difference in incidence between two groups
"""
# Method for single survey
function incprops(prev::Float64,
    σ_prev::Float64,
    prevR::Float64,
    σ_prevR::Float64,
    mdri::Float64,
    σ_mdri::Float64,
    frr::Float64,
    σ_frr::Float64;
    covar::Float64 = 0.0, # covariance of prev and prevR
    T = 730.5, # in same units as MDRI
    timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
    bs::Int64 = 0,
    gibbs = false,
    bs_numbers = false,
    bs_numbers_n::AbstractVector{Int64} = [0, 0],
    α::Float64 = 0.05,
    per::Int64 = 1)

    # convert to estimation unit
    mdri = mdri / timeconversion
    σ_mdri = σ_mdri / timeconversion
    T = T / timeconversion

    pe = kassanjee(prev, prevR, mdri, frr, T) * per

    if bs == 0 && σ_prev == 0.0
        @warn "σ_prev of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if bs == 0 && σ_prevR == 0.0
        @warn "σ_prevR of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if bs == 0 && σ_mdri == 0.0
        @warn "σ_mdri of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if bs == 0 && σ_frr == 0.0
        @warn "σ_frr of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if bs > 0 && σ_prev == 0.0
        @warn "σ_prev of zero supplied. Set to 0.0000000001."
        σ_prev = 0.0000000001
    end
    if bs > 0 && σ_prevR == 0.0
        @warn "σ_prevR of zero supplied. Set to 0.0000000001."
        σ_prevR = 0.0000000001
    end
    if bs > 0 && σ_mdri == 0.0
        @warn "σ_mdri of zero supplied. Set to 0.0000000001."
        σ_mdri = 0.0000000001
    end
    if bs > 0 && σ_frr == 0.0
        @warn "σ_frr of zero supplied. Set to 0.0000000001."
        σ_frr = 0.0000000001
    end

    if bs == 0 && bs_numbers
        @error "Cannot bootstrap numbers if bootstrapping is not being performed"
    end

    if bs_numbers && any(x->x==0, bs_numbers_n)
        @error "Cannot bootstrap numbers if number of trials is zero"
    end

    if bs_numbers && bs_numbers_n[1] == 0
        @warn "Set n1 to 1"
        bs_numbers_n[1] = 1
    end

    if bs_numbers && bs_numbers_n[2] == 0
        @warn "Set n2 to 1"
        bs_numbers_n[2] = 1
    end

    if bs == 0
        σ, σ_infSS = σ_dm(prev, prevR, mdri, frr, T, σ_prev, σ_prevR, σ_mdri, σ_frr) .* per
        ci = Distributions.quantile.(Distributions.Normal(pe, σ), [α/2, 1-α/2]) # max.(Distributions.quantile.(Distributions.Normal(pe, σ), [α/2, 1-α/2]),0)
        return (I = pe, CI = ci, σ = σ, RSE = σ/abs(pe))

    # Manual implementation of truncated normal distribution
    elseif bs > 0 && !bs_numbers && gibbs
        if covar < 0.0
            @warn "Covariance of prev and prevR cannot be negative, set to 0.0"
            covar = 0.0
        elseif covar == 0.0
            r = hcat(
                    rtnorm_gibbs(bs, prev, σ_prev, 0.0, 1.0),
                    rtnorm_gibbs(bs, prevR, σ_prevR, 0.0, 1.0),
                    rtnorm_gibbs(bs, mdri, σ_mdri, 0.0, Inf),
                    rtnorm_gibbs(bs, frr, σ_frr, 0.0, 1.0)
                    )
        elseif covar > 0.0
            # Use rejection sampling for correlated prev and prevR
            µ = [prev, prevR, mdri, frr]
            Σ = [σ_prev^2 covar 0 0 ; covar σ_prevR^2 0 0 ; 0 0 σ_mdri^2 0 ; 0 0 0 σ_frr^2]
            r = rtmvnorm(bs, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, Inf, 1.0])
        end

        bs_incidence = kassanjee.(r[:,1], r[:,2], r[:,3], r[:,4], T) .* per
        σ = Statistics.std(bs_incidence)
        ci = Statistics.quantile(bs_incidence, [α/2, 1-α/2]) # max.(Statistics.quantile(bs_incidence, [α/2, 1-α/2]),0)
        cov_prev_I = Statistics.cov([r[:,1] bs_incidence])
        cor_prev_I = Statistics.cor([r[:,1] bs_incidence])
        return (I = pe, CI = ci, σ = σ, RSE = σ/abs(pe), cov_prev_I = cov_prev_I, cor_prev_I = cor_prev_I)

    # Use Distributions.jl truncated normal
    elseif bs > 0 && !bs_numbers && !gibbs
        if covar < 0.0
            @warn "Covariance of prev and prevR cannot be negative, set to 0.0"
            covar = 0.0
        elseif covar == 0.0
            r = hcat(
                    rand(Distributions.truncated(Distributions.Normal(prev, σ_prev), 0.0, 1.0), bs),
                    rand(Distributions.truncated(Distributions.Normal(prevR, σ_prevR), 0.0, 1.0), bs),
                    rand(Distributions.truncated(Distributions.Normal(mdri, σ_mdri), 0.0, Inf), bs),
                    rand(Distributions.truncated(Distributions.Normal(frr, σ_frr), 0.0, 1.0), bs)
                    )
        elseif covar > 0.0
            # Use rejection sampling for correlated prev and prevR
            µ = [prev, prevR, mdri, frr]
            Σ = [σ_prev^2 covar 0 0 ; covar σ_prevR^2 0 0 ; 0 0 σ_mdri^2 0 ; 0 0 0 σ_frr^2]
            r = rtmvnorm(bs, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, Inf, 1.0])
        end

        bs_incidence = kassanjee.(r[:,1], r[:,2], r[:,3], r[:,4], T) .* per
        σ = Statistics.std(bs_incidence)
        ci = Statistics.quantile(bs_incidence, [α/2, 1-α/2]) # max.(Statistics.quantile(bs_incidence, [α/2, 1-α/2]),0)
        cov_prev_I = Statistics.cov([r[:,1] bs_incidence])
        cor_prev_I = Statistics.cor([r[:,1] bs_incidence])
        return (I = pe, CI = ci, σ = σ, RSE = σ/abs(pe), cov_prev_I = cov_prev_I, cor_prev_I = cor_prev_I)

    elseif bs > 0 && bs_numbers
        @warn "Covariance between prevalence and prevalence of recency assumed zero"
        dprev = Distributions.Binomial(bs_numbers_n[1], prev)
        dprevR = Distributions.Binomial(bs_numbers_n[2], prevR)
        dmdri = Distributions.truncated(Distributions.Normal(mdri, σ_mdri), 0, Inf)
        dfrr = Distributions.truncated(Distributions.Normal(frr, σ_frr), 0, 1)
        npos = rand(dprev, bs)
        prevs =  npos ./ bs_numbers_n[1]
        nr = rand(dprevR, bs)
        prevRs = nr ./ max.(min.(bs_numbers_n[2], npos), 1)
        mdris = rand(dmdri, bs)
        frrs = rand(dfrr, bs)
        bs_incidence = kassanjee.(prevs, prevRs, mdris, frrs, T) .* per
        σ = Statistics.std(bs_incidence)
        ci = Statistics.quantile(bs_incidence, [α/2, 1-α/2]) # max.(Statistics.quantile(bs_incidence, [α/2, 1-α/2]),0)
        cov_prev_I = Statistics.cov([prevs bs_incidence])
        cor_prev_I = Statistics.cor([prevs bs_incidence])
        return (I = pe, CI = ci, σ = σ, RSE = σ/abs(pe), cov_prev_I = cov_prev_I, cor_prev_I = cor_prev_I)
    end
end

"""
    incdif(prev, σ_prev, prevR, σ_prevR, mdri, σ_mdri, frr, σ_frr; covar, kwargs...)

Estimate the difference in HIV incidence between two groups or time points.

This function calculates the difference between two incidence estimates (Group 1 - Group 2)
and provides statistical inference including confidence intervals and hypothesis testing.
Uses the Kassanjee incidence estimator for each group.

# Arguments
## Required Parameters (Vectors for Groups 1 and 2)
- `prev::Vector{Float64}`: Prevalences for [group1, group2]
- `σ_prev::Vector{Float64}`: Standard errors of prevalences
- `prevR::Vector{Float64}`: Recency prevalences for [group1, group2]
- `σ_prevR::Vector{Float64}`: Standard errors of recency prevalences
- `mdri::Float64`: Mean duration of recent infection (same for both groups)
- `σ_mdri::Float64`: Standard error of MDRI
- `frr::Float64`: False-recent rate (same for both groups)
- `σ_frr::Float64`: Standard error of FRR

## Required Keyword Arguments
- `covar::Vector{Float64}`: Covariances [cov(prev1,prevR1), cov(prev2,prevR2)]

## Optional Keyword Arguments
- `T::Float64=730.5`: Recency test time cutoff (typically in days)
- `timeconversion::Float64=365.25`: Convert MDRI/T units to incidence time units
- `bs::Int64=0`: Number of bootstrap iterations (0 = Delta method only)
- `output_bs::Bool=false`: Return bootstrap samples for further analysis
- `bs_numbers::Bool=false`: Bootstrap from counts rather than proportions
- `bs_numbers_n::Vector{Int64}=[0,0,0,0]`: Sample sizes [n1_total, n1_pos, n2_total, n2_pos]
- `α::Float64=0.05`: Significance level (adjusted by Bonferroni correction if specified)
- `bonf_cor::Int64=1`: Bonferroni correction factor for multiple comparisons
- `per::Int64=1`: Multiplier for incidence rates

# Returns
A named tuple containing:
- `Δ::Float64`: Point estimate of incidence difference (Group 1 - Group 2)
- `CI::Vector{Float64}`: Confidence interval [lower, upper]
- `σ::Float64`: Standard error of the difference
- `RSE::Float64`: Relative standard error (σ/|Δ|)
- `p::Float64`: Two-sided p-value for H₀: Δ = 0
- `bs_difs::Vector{Float64}`: Bootstrap differences (only if `output_bs=true`)

# Statistical Methods
1. **Delta Method** (`bs=0`): Analytical variance using error propagation
2. **Bootstrap** (`bs>0`): Resampling with correlation structure preserved
3. **Hypothesis Test**: Normal approximation Z-test for difference

# Notes
- Negative incidence estimates are set to 0.0 before computing differences
- Bonferroni correction adjusts α for multiple comparisons (α_adj = α / bonf_cor)
- When `bs_numbers=true`, bootstraps from binomial distributions for finite samples
- Bootstrap samples preserve correlation within groups but assume independence between groups

# Examples
```julia
# Example 1: Compare incidence between two regions
result = incdif(
    [0.20, 0.15],      # Prevalences: Region 1 = 20%, Region 2 = 15%
    [0.015, 0.012],    # Standard errors
    [0.10, 0.08],      # Recency prevalences
    [0.02, 0.015],     # Standard errors
    130.0, 15.0,       # MDRI = 130 ± 15 days
    0.01, 0.005,       # FRR = 1% ± 0.5%
    covar = [0.0, 0.0] # No covariance within groups
)
println("Incidence difference: \$(result.Δ)")
println("95% CI: \$(result.CI)")
println("p-value: \$(result.p)")

# Example 2: Multiple comparisons with Bonferroni correction
result_bonf = incdif(
    [0.20, 0.15], [0.015, 0.012],
    [0.10, 0.08], [0.02, 0.015],
    130.0, 15.0, 0.01, 0.005,
    covar = [0.0, 0.0],
    bonf_cor = 3  # Comparing 3 pairs, adjust α
)

# Example 3: Bootstrap with output for custom analysis
result_bs = incdif(
    [0.20, 0.15], [0.015, 0.012],
    [0.10, 0.08], [0.02, 0.015],
    130.0, 15.0, 0.01, 0.005,
    covar = [0.0, 0.0],
    bs = 10000,
    output_bs = true
)
# Access bootstrap samples
histogram(result_bs.bs_difs)
```

# References
Kassanjee R, et al. (2012). A new general biomarker-based incidence estimator.
Epidemiology, 23(5), 721-728.

# See Also
- `incprops`: Estimate single incidence
- `inccounts`: Estimate incidence from counts
"""
# Incidence difference
function incdif(prev::AbstractVector{Float64},
    σ_prev::AbstractVector{Float64},
    prevR::AbstractVector{Float64},
    σ_prevR::AbstractVector{Float64},
    mdri::Float64,
    σ_mdri::Float64,
    frr::Float64,
    σ_frr::Float64;
    covar::AbstractVector{Float64}, # covariances of prev and prevR
    T = 730.5, # in same units as MDRI
    timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
    bs::Int64 = 0,
    output_bs = false,
    bs_numbers = false,
    bs_numbers_n::AbstractVector{Int64} = [0, 0, 0, 0],
    α::Float64 = 0.05,
    bonf_cor::Int64 = 1,
    per::Int64 = 1)

    if bonf_cor < 1
        @error "Bonferroni correction only possible with positive integers for number of comparisons"
    end

    if bonf_cor > 1
        α = α / bonf_cor
    end

    if bs == 0 && σ_mdri == 0
        @warn "σ_mdri of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if bs == 0 && σ_frr == 0
        @warn "σ_frr of zero supplied. Variance of incidence estimate likely incorrect."
    end

    if bs > 0 && σ_mdri == 0
        @warn "σ_mdri of zero supplied. Set to 0.0000000001."
        σ_mdri = 0.0000000001
    end
    if bs > 0 && σ_frr == 0
        @warn "σ_frr of zero supplied. Set to 0.0000000001."
        σ_frr = 0.0000000001
    end


    if bs == 0 && σ_prev[1] == 0
        @warn "σ_prev of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if bs == 0 && σ_prevR[1] == 0
        @warn "σ_prevR of zero supplied. Variance of incidence estimate likely incorrect."
    end

    if bs == 0 && σ_prev[2] == 0
        @warn "σ_prev of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if bs == 0 && σ_prevR[2] == 0
        @warn "σ_prevR of zero supplied. Variance of incidence estimate likely incorrect."
    end

    if bs > 0 && σ_prev[1] == 0
        @warn "σ_prev of zero supplied. Set to 0.0000000001."
        σ_prev[1] = 0.0000000001
    end
    if bs > 0 && σ_prevR[1] == 0
        @warn "σ_prevR of zero supplied. Set to 0.0000000001."
        σ_prevR[1] = 0.0000000001
    end

    if bs > 0 && σ_prev[2] == 0
        @warn "σ_prev of zero supplied. Set to 0.0000000001."
        σ_prev[2] = 0.0000000001
    end
    if bs > 0 && σ_prevR[2] == 0
        @warn "σ_prevR of zero supplied. Set to 0.0000000001."
        σ_prevR[2] = 0.0000000001
    end

    # convert to estimation unit
    mdri = mdri / timeconversion
    σ_mdri = σ_mdri / timeconversion
    T = T / timeconversion

    I_1 = kassanjee(prev[1], prevR[1], mdri, frr, T) * per
    if I_1 < 0.0
        @warn "Negative point estimate set to 0.0 for Δ calculation"
        I_1 = 0.0
    end
    I_2 = kassanjee(prev[2], prevR[2], mdri, frr, T) * per
    if I_2 < 0.0
        @warn "Negative point estimate set to 0.0 for Δ calculation"
        I_2 = 0.0
    end
    pe = I_1 - I_2

    if σ_prev == 0
        @warn "σ_prev of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if σ_prevR == 0
        @warn "σ_prevR of zero supplied. Variance of incidence estimate likely incorrect."
    end
    if σ_mdri == 0
        @warn "σ_mdri of zero supplied."
    end
    if σ_frr == 0
        @warn "σ_frr of zero supplied."
    end

    if bs == 0 && bs_numbers
        @error "Cannot bootstrap numbers if bootstrapping is not being performed"
    end

    if bs_numbers && any(x->x==0, bs_numbers_n)
        @error "Cannot bootstrap numbers if number of trials is zero"
    end

    if bs_numbers && bs_numbers_n[1] == 0
        @warn "Set n1 to 1"
        bs_numbers_n[1] = 1
    end

    if bs_numbers && bs_numbers_n[2] == 0
        @warn "Set n2 to 1"
        bs_numbers_n[2] = 1
    end

    if bs_numbers && length(bs_numbers_n) == 4 && bs_numbers_n[3] == 0
        @warn "Set n3 to 1"
        bs_numbers_n[3] = 1
    end

    if bs_numbers && length(bs_numbers_n) == 4 && bs_numbers_n[4] == 0
        @warn "Set n4 to 1"
        bs_numbers_n[4] = 1
    end

    if bs == 0
        σ = σ_Δ_dm(prev, prevR, mdri, frr, T, σ_prev, σ_prevR, σ_mdri, σ_frr) * per
        ci = Distributions.quantile.(Distributions.Normal(pe, σ), [α/2, 1-α/2])
        p = Distributions.cdf(Distributions.Normal(), -abs(pe)/σ) * 2
        return (Δ = pe, CI = ci, σ = σ, RSE = σ/abs(pe), p = p)

    elseif bs > 0 && !bs_numbers
        if covar[1] < 0.0
            @warn "Covariance of prev and prevR cannot be negative, set to 0.0"
            covar[1] = 0.0
        end
        if covar[2] < 0.0
            @warn "Covariance of prev and prevR cannot be negative, set to 0.0"
            covar[2] = 0.0
        end
        if covar[1] == 0.0 && covar[2] == 0.0
            r = hcat(
                    rand(Distributions.truncated(Distributions.Normal(prev[1], σ_prev[1]), 0.0, 1.0), bs),
                    rand(Distributions.truncated(Distributions.Normal(prevR[1], σ_prevR[1]), 0.0, 1.0), bs),
                    rand(Distributions.truncated(Distributions.Normal(prev[2], σ_prev[2]), 0.0, 1.0), bs),
                    rand(Distributions.truncated(Distributions.Normal(prevR[2], σ_prevR[2]), 0.0, 1.0), bs),
                    rand(Distributions.truncated(Distributions.Normal(mdri, σ_mdri), 0.0, Inf), bs),
                    rand(Distributions.truncated(Distributions.Normal(frr, σ_frr), 0.0, 1.0), bs)
                    )
        elseif covar[1] > 0.0 || covar[2] > 0.0
            @error "Truncated 6D multivariate normal with covariance not supported. Rejection sampling is limited to 4 dimensions. Set covar=[0.0, 0.0] to use independent Gibbs sampling."
            µ = [prev[1], prevR[1], prev[2], prevR[2], mdri, frr]
            Σ = [σ_prev[1]^2 covar[1] 0 0 0 0; covar[1] σ_prevR[1]^2 0 0 0 0 ; 0 0 σ_prev[2]^2 covar[2] 0 0 ; 0 0 covar[2] σ_prevR[2]^2 0 0 ; 0 0 0 0 σ_mdri^2 0 ; 0 0 0 0 0 σ_frr^2]
            r = rtmvnorm(bs, µ, Σ, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0, Inf, 1.0])
        end

        bs_incidence_1 = kassanjee.(r[:,1], r[:,2], r[:,5], r[:,6], T) .* per
        bs_incidence_2 = kassanjee.(r[:,3], r[:,4], r[:,5], r[:,6], T) .* per
        if any(x-> x < 0, bs_incidence_1) || any(x-> x < 0, bs_incidence_2)
            @warn "Negative bootstrapped incidence estimates set to 0.0 for difference calculations"
        end
        bs_difs = max.(bs_incidence_1,0.0) .- max.(bs_incidence_2,0.0)
        σ = Statistics.std(bs_difs)
        ci = Statistics.quantile(bs_difs, [α/2, 1-α/2])
        p = Distributions.cdf(Distributions.Normal(), -abs(pe)/σ) * 2
        # These numbers differ from R - bug?
        #cov_prev_I = [Statistics.cov([r[:,1] bs_incidence_1]) Statistics.cov([r[:,3] bs_incidence_2])]
        #cor_prev_I = [Statistics.cor([r[:,1] bs_incidence_1]) Statistics.cor([r[:,3] bs_incidence_2])]
        if output_bs
            return (Δ = pe, CI = ci, σ = σ, RSE = σ/abs(pe), p = p, bs_difs = bs_difs) #, cov_prev_I = cov_prev_I, cor_prev_I = cor_prev_I #ptest = ptest,
        else
            return (Δ = pe, CI = ci, σ = σ, RSE = σ/abs(pe), p = p)
        end

    elseif bs > 0 && bs_numbers
        @warn "Covariance between prevalence and prevalence of recency assumed zero"
        dprev_1 = Distributions.Binomial(bs_numbers_n[1], prev[1])
        dprevR_1 = Distributions.Binomial(bs_numbers_n[2], prevR[1])
        dprev_2 = Distributions.Binomial(bs_numbers_n[3], prev[2])
        dprevR_2 = Distributions.Binomial(bs_numbers_n[4], prevR[2])
        dmdri = Distributions.truncated(Distributions.Normal(mdri, σ_mdri), 0, Inf)
        dfrr = Distributions.truncated(Distributions.Normal(frr, σ_frr), 0, 1)
        npos_1 = rand(dprev_1, bs)
        prevs_1 =  npos_1 ./ bs_numbers_n[1]
        nr_1 = rand(dprevR_1, bs)
        prevRs_1 = nr_1 ./ bs_numbers_n[2]
        npos_2 = rand(dprev_2, bs)
        prevs_2 =  npos_2 ./ bs_numbers_n[3]
        nr_2 = rand(dprevR_2, bs)
        prevRs_2 = nr_2 ./ bs_numbers_n[4]
        mdris = rand(dmdri, bs)
        frrs = rand(dfrr, bs)
        bs_incidence_1 = kassanjee.(prevs_1, prevRs_1, mdris, frrs, T) .* per
        bs_incidence_2 = kassanjee.(prevs_2, prevRs_2, mdris, frrs, T) .* per
        if any(x->x<0, bs_incidence_1) || any(x->x<0, bs_incidence_2)
            @warn "Negative bootstrapped incidence estimates set to 0.0 for difference calculations"
        end
        bs_difs = max.(bs_incidence_1,0.0) .- max.(bs_incidence_2,0.0)
        σ = Statistics.std(bs_difs)
        ci = Statistics.quantile(bs_difs, [α/2, 1-α/2])
        bs_difs_abs = abs.(bs_difs)
        p = Distributions.cdf(Distributions.Normal(), -abs(pe)/σ) * 2
        if output_bs
            return (Δ = pe, CI = ci, σ = σ, RSE = σ/abs(pe), p = p, bs_difs = bs_difs) #, cov_prev_I = cov_prev_I, cor_prev_I = cor_prev_I #ptest = ptest,
        else
            return (Δ = pe, CI = ci, σ = σ, RSE = σ/abs(pe), p = p)
        end
    end
end

# Method for multiple surveys using one test
# function incprops(prev::AbstractVector{Float64},
#     prevR::AbstractVector{Float64},
#     mdri::Float64, # in days
#     frr::Float64; # variance-covariance matrix for prev and prevR
#     σ_mdri::Float64 = 0.0,
#     σ_frr::Float64 = 0.0,
#     covar::Array{Float64,2} = Matrix{Float64}(I, size(prev)[1], size(prev)[1]),
#     T = 730.5, # in same units as MDRI
#     timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
#     bs::Int64 = 0,
#     α::Float64 = 0.05)
#
#     # convert to estimation unit
#     mdri = mdri / timeconversion
#     σ_mdri = σ_mdri / timeconversion
#     T = T / timeconversion
#
#     pe = kassanjee.(prev, prevR, mdri, frr, T)
#
#     for i in 1:size(prev)[1]
#         if bs == 0 && σ_prev > 0 && σ_prevR > 0 && σ_mdri > 0 && σ_frr > 0
#     end
#     if bs == 0 && σ_prev > 0 && σ_prevR > 0 && σ_mdri > 0 && σ_frr > 0
#         σ = σ_dm(prev, prevR, mdri, frr, T, σ_prev, σ_prevR, σ_mdri, σ_frr) * per
#         ci = max.(Distributions.quantile.(Distributions.Normal(pe, σ), [α/2, 1-α/2]),0)
#         return pe, σ, ci
#     else
#         return incidence_pe
#     end
# end

# Method for multiple surveys using multiple tests
# function incprops(prev::AbstractVector{Float64},
#     prevR::AbstractVector{Float64},
#     mdri::AbstractVector{Float64}, # in days
#     frr::AbstractVector{Float64};
#     σ_mdri::AbstractVector{Float64} = repeat([0.0],size(mdri)[1]),
#     σ_frr::AbstractVector{Float64} = repeat([0.0],size(frr)[1]),
#     covar::Array{Float64,2} = Matrix{Float64}(I, size(prev)[1], size(prev)[1]),
#     T::AbstractVector{Float64} = repeat([730.5],size(mdri)[1]), # in same unit as MDRI
#     timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
#     bs::Int64 = 0,
#     α::Float64 = 0.05)
#
#     # convert to estimation unit
#     mdri = mdri / timeconversion
#     σ_mdri = σ_mdri / timeconversion
#     T = T / timeconversion
#
#     incidence_pe = kassanjee.(prev, prevR, mdri, frr, T)
#     return incidence_pe
# end

"""
    inccounts(n, npos, ntestR, nR, mdri, frr; kwargs...)
    inccounts(n::Vector, npos::Vector, ntestR::Vector, nR::Vector, mdri, frr; kwargs...)

Estimate HIV incidence from survey count data (single or multiple surveys).

This is a convenience wrapper around `incprops` that first estimates prevalences from counts
using the `prevalence` function, then calculates incidence. Automatically handles standard
error propagation from count data.

# Arguments (Single Survey Method)
## Required Parameters
- `n::Int64`: Total sample size tested for HIV
- `npos::Int64`: Number testing HIV-positive
- `ntestR::Int64`: Number of HIV-positive individuals tested for recency
- `nR::Int64`: Number testing recent among those tested
- `mdri::Float64`: Mean duration of recent infection (typically in days)
- `frr::Float64`: False-recent rate (proportion)

## Keyword Arguments
- `de_npos::Float64=1.0`: Design effect for HIV prevalence estimate
- `de_nR::Float64=1.0`: Design effect for recency prevalence estimate
- `σ_mdri::Float64=0.0`: Standard error of MDRI
- `σ_frr::Float64=0.0`: Standard error of FRR
- `cov::Float64=0.0`: Covariance between prevalence and recency prevalence
- `T::Float64=730.5`: Recency test time cutoff
- `timeconversion::Float64=365.25`: Time unit conversion factor
- `bs::Int64=0`: Number of bootstrap iterations
- `α::Float64=0.05`: Significance level for confidence intervals
- `per::Int64=1`: Multiplier for incidence rate

# Arguments (Multiple Surveys Method)
Same as single survey but with vectors for survey-specific counts:
- `n::Vector{Int64}`: Sample sizes for each survey
- `npos::Vector{Int64}`: HIV-positive counts
- `ntestR::Vector{Int64}`: Tested for recency counts
- `nR::Vector{Int64}`: Recent infection counts
- `de_npos::Vector{Float64}`: Design effects for prevalence
- `de_nR::Vector{Float64}`: Design effects for recency
- `cov::Matrix{Float64}`: Variance-covariance matrix

# Returns
Same as `incprops`: Named tuple with I, CI, σ, RSE, and optionally cov_prev_I, cor_prev_I

# Workflow
1. Calculate prevalence: `prev = npos / n` with standard error
2. Calculate recency prevalence: `prevR = nR / ntestR` with standard error
3. Call `incprops` with calculated proportions and standard errors

# Examples
```julia
# Example 1: Single survey from counts
result = inccounts(
    5000,      # 5,000 people tested for HIV
    1000,      # 1,000 HIV-positive
    900,       # 900 of those tested for recency
    90,        # 90 tested recent
    130.0,     # MDRI = 130 days
    0.01,      # FRR = 1%
    σ_mdri = 15.0,
    σ_frr = 0.005,
    bs = 10000
)
println("Incidence: \$(result.I) per year")

# Example 2: With design effects for complex survey
result_de = inccounts(
    5000, 1000, 900, 90,
    130.0, 0.01,
    de_npos = 1.5,    # Clustering effect
    de_nR = 1.3,
    σ_mdri = 15.0,
    σ_frr = 0.005
)

# Example 3: Multiple surveys (returns vector of results)
results = inccounts(
    [5000, 4500, 5200],  # Three surveys
    [1000, 900, 1100],   # HIV-positive counts
    [900, 850, 1000],    # Tested for recency
    [90, 85, 105],       # Recent infections
    130.0, 0.01,
    de_npos = [1.0, 1.0, 1.0],
    de_nR = [1.0, 1.0, 1.0],
    σ_mdri = 15.0,
    σ_frr = 0.005
)
```

# Notes
- Uses normal approximation for prevalence standard errors
- Automatically accounts for finite population sampling
- Design effects inflate standard errors to account for clustering
- For `ntestR < npos`, only a subset was tested for recency
- Warnings issued if counts are too small for normal approximation

# See Also
- `prevalence`: Calculate prevalence from counts
- `incprops`: Estimate incidence from proportions (underlying function)
- `incdif`: Compare incidence between groups
"""
# Single survey
function inccounts(n::Int64,
    npos::Int64,
    ntestR::Int64,
    nR::Int64,
    mdri::Float64, # in days
    frr::Float64;
    de_npos::Float64 = 1.0,
    de_nR::Float64 = 1.0,
    σ_mdri::Float64 = 0.0,
    σ_frr::Float64 = 0.0,
    covar::Union{Float64, Nothing} = nothing,  # covariance of prev and prevR
    cov::Union{Float64, Nothing} = nothing,    # DEPRECATED: use covar instead
    T = 730.5, # in same units as MDRI
    timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
    bs::Int64 = 0,
    α::Float64 = 0.05,
    per::Int64 = 1)

    # Handle deprecated cov parameter
    if cov !== nothing && covar !== nothing
        if cov != covar
            error("Both 'cov' and 'covar' were provided with different values. Please use only 'covar' (cov is deprecated).")
        end
        @warn "The 'cov' argument is deprecated and will be removed in the next version. Please use 'covar' instead."
        covar_value = covar
    elseif cov !== nothing
        @warn "The 'cov' argument is deprecated and will be removed in the next version. Please use 'covar' instead."
        covar_value = cov
    elseif covar !== nothing
        covar_value = covar
    else
        covar_value = 0.0  # Default value
    end

    # compute prevalences
    prev, σ_prev  = prevalence(npos, n, de_npos)
    prevR, σ_prevR = prevalence(nR, ntestR, de_nR)

    return incprops(prev, σ_prev, prevR, σ_prevR, mdri, σ_mdri, frr, σ_frr,
                            covar = covar_value, T = T, timeconversion = timeconversion,
                            bs = bs, α = α, per = per)
end

# Multiple surveys - single test
function inccounts(n::AbstractVector{Int64},
    npos::AbstractVector{Int64},
    ntestR::AbstractVector{Int64},
    nR::AbstractVector{Int64},
    mdri::Float64, # in days
    frr::Float64;
    de_npos::AbstractVector{Float64} = repeat([1.0],size(n)[1]),
    de_nR::AbstractVector{Float64} = repeat([1.0],size(n)[1]),
    σ_mdri::Float64 = 0.0,
    σ_frr::Float64 = 0.0,
    covar::Union{Array{Float64,2}, Nothing} = nothing,  # covariance matrix of prev and prevR
    cov::Union{Array{Float64,2}, Nothing} = nothing,    # DEPRECATED: use covar instead
    T = 730.5, # in same units as MDRI
    timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
    bs::Int64 = 0,
    α::Float64 = 0.05,
    per::Int64 = 1)

    # compute prevalences
    prev, σ_prev  = prevalence(npos, n, de_npos)
    prevR, σ_prevR = prevalence(nR, ntestR, de_nR)

    # Handle deprecated cov parameter
    if cov !== nothing && covar !== nothing
        if cov != covar
            error("Both 'cov' and 'covar' were provided with different values. Please use only 'covar' (cov is deprecated).")
        end
        @warn "The 'cov' argument is deprecated and will be removed in the next version. Please use 'covar' instead."
        covar_value = covar
    elseif cov !== nothing
        @warn "The 'cov' argument is deprecated and will be removed in the next version. Please use 'covar' instead."
        covar_value = cov
    elseif covar !== nothing
        covar_value = covar
    else
        # Default: identity matrix (no correlation between surveys)
        covar_value = Matrix{Float64}(I, size(prev)[1], size(prev)[1])
    end

    return incprops(prev, σ_prev, prevR, σ_prevR, mdri, σ_mdri, frr, σ_frr,
                            covar = covar_value, T = T, timeconversion = timeconversion,
                            bs = bs, α = α, per = per)
end

end
