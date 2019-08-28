module Incidence

export prevalence, rtmvnorm, incprops, inccounts

import LinearAlgebra.I
import Distributions
import Statistics

function prevalence(pos, n, de = 1) #, f = 1
    p = pos/n
    σ = sqrt( (p * (1 - p)) / n ) * de #* sqrt(1 - f)
    return p, σ
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

function rtmvnorm(n::Int64,
    µ::AbstractVector{Float64},
    Σ::Array{Float64,2},
    lower::AbstractVector{Float64} = [-Inf, -Inf, -Inf, -Inf],
    upper::AbstractVector{Float64} = [Inf, Inf, Inf, Inf])

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
    α::Float64 = 0.05,
    per::Int64 = 1)

    # convert to estimation unit
    mdri = mdri / timeconversion
    σ_mdri = σ_mdri / timeconversion
    T = T / timeconversion

    pe = kassanjee(prev, prevR, mdri, frr, T) * per

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

    if bs == 0
        σ, σ_infSS = σ_dm(prev, prevR, mdri, frr, T, σ_prev, σ_prevR, σ_mdri, σ_frr) .* per
        ci = max.(Distributions.quantile.(Distributions.Normal(pe, σ), [α/2, 1-α/2]),0)
        return (I = pe, CI = ci, σ = σ, RSE = σ/pe)

    # Manual implementation of truncated normal distribution
    elseif bs > 0 #&& covar > 0.0
        µ = [prev, prevR, mdri, frr]
        Σ = [σ_prev^2 covar 0 0 ; covar σ_prevR^2 0 0 ; 0 0 σ_mdri^2 0 ; 0 0 0 σ_frr^2]
        r = rtmvnorm(bs, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, Inf, 1.0])
        bs_incidence = kassanjee.(r[:,1], r[:,2], r[:,3], r[:,4], T) * per
        σ = Statistics.std(bs_incidence)
        ci = max.(Statistics.quantile(bs_incidence, [α/2, 1-α/2]),0)
        # These numbers differ from R - bug?
        cov_prev_I = Statistics.cov([r[:,1] bs_incidence])
        cor_prev_I = Statistics.cor([r[:,1] bs_incidence])
        return (I = pe, CI = ci, σ = σ, RSE = σ/pe, cov_prev_I = cov_prev_I, cor_prev_I = cor_prev_I)
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
function incprops(prev::AbstractVector{Float64},
    prevR::AbstractVector{Float64},
    mdri::AbstractVector{Float64}, # in days
    frr::AbstractVector{Float64};
    σ_mdri::AbstractVector{Float64} = repeat([0.0],size(mdri)[1]),
    σ_frr::AbstractVector{Float64} = repeat([0.0],size(frr)[1]),
    covar::Array{Float64,2} = Matrix{Float64}(I, size(prev)[1], size(prev)[1]),
    T::AbstractVector{Float64} = repeat([730.5],size(mdri)[1]), # in same unit as MDRI
    timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
    bs::Int64 = 0,
    α::Float64 = 0.05)

    # convert to estimation unit
    mdri = mdri / timeconversion
    σ_mdri = σ_mdri / timeconversion
    T = T / timeconversion

    incidence_pe = kassanjee.(prev, prevR, mdri, frr, T)
    return incidence_pe
end

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
    cov::Float64 = 0.0, # covariance of prev and prevR
    T = 730.5, # in same units as MDRI
    timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
    bs::Int64 = 0,
    α::Float64 = 0.05,
    per::Int64 = 1)

    # compute prevalences
    prev, σ_prev  = prevalence(npos, n, de_npos)
    prevR, σ_prevR = prevalence(nR, ntestR, de_nR)

    return incprops(prev, prevR, mdri, frr, σ_prev = σ_prev,
                            σ_prevR = σ_prevR, σ_mdri = σ_mdri, σ_frr = σ_frr,
                            cov = 0.0, T = T, timeconversion = timeconversion,
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
    cov::Array{Float64,2} = Matrix{Float64}(I, size(prev)[1], size(prev)[1]), # covariance of prev and prevR
    T = 730.5, # in same units as MDRI
    timeconversion = 365.25, # to convert from unit in which MDRI and T is specified to unit of incidence
    bs::Int64 = 0,
    α::Float64 = 0.05,
    per::Int64 = 1)

    # compute prevalences
    prev, σ_prev  = prevalence(npos, n, de_npos)
    prevR, σ_prevR = prevalence(nR, ntestR, de_nR)

    return incprops(prev, prevR, mdri, frr, σ_prev = σ_prev,
                            σ_prevR = σ_prevR, σ_mdri = σ_mdri, σ_frr = σ_frr,
                            cov = 0.0, T = T, timeconversion = timeconversion,
                            bs = bs, α = α, per = per)
end



end
