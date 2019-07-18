module Incidence

export prevalence, incprops, inccounts

import LinearAlgebra.I
using Distributions
using Statistics

import Distributions.Normal
import Distributions.Truncated
import Distributions.rand

function prevalence(pos, n, de = 1) #, f = 1
    p = pos/n
    σ = sqrt( (p * (1 - p)) / n ) * de #* sqrt(1 - f)
    return p, σ
end

function kassanjee(prev::Float64, prevR::Float64, mdri::Float64, frr::Float64, T)
    return (prev * (prevR - frr)) / ((1 - prev) * (mdri - frr * T))
end

# Method for single survey
function incprops(prev::Float64,
    prevR::Float64,
    mdri::Float64,
    frr::Float64;
    σ_prev::Float64 = 0.0,
    σ_prevR::Float64 = 0.0,
    σ_mdri::Float64 = 0.0,
    σ_frr::Float64 = 0.0,
    cov::Float64 = 0.0, # covariance of prev and prevR
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

    if bs == 0
        fot_prev = (prevR - frr)/(((1 - prev)^2) * (mdri - frr * T))  #E.G. d(I)/d(P_H)
        fot_prevR = prev/((1 - prev) * (mdri - frr * T))
        fot_mdri = (frr * prev - prevR * prev)/((1 - prev) * ((mdri - frr * T)^2))
        fot_frr = (prev * (T * prev - mdri))/((1 - prev) * ((mdri - frr * T)^2))

        dm_sd = sqrt(fot_prev^2 * σ_prev^2 + fot_prevR^2 *
            σ_prevR^2 + fot_mdri^2 * σ_mdri^2 + fot_frr^2 * σ_frr^2) * per

        dm_ci = quantile.(Normal(pe, dm_sd), [α/2, 1-α/2])

        return pe, dm_sd, dm_ci

    elseif bs > 0
        print("Warning: In this implementation, prev and prevR are assumed independent. Variance-covariance matrix is:\n")
        vcovmat = [σ_prev^2 0 0 0 ; 0 σ_prevR^2 0 0 ; 0 0 σ_mdri^2 0 ; 0 0 0 σ_frr^2]
        display(vcovmat)
        if σ_prev == 0
            print("Warning: bootstrapping with σ_prev = 0")
        end
        if σ_prevR == 0
            print("Warning: bootstrapping with σ_prevR = 0")
        end
        if σ_mdri == 0
            print("Warning: bootstrapping with σ_mdri = 0")
        end
        if σ_frr == 0
            print("Warning: bootstrapping with σ_frr = 0")
        end
        # let's fist do it with independent prev and prevR
        tn_prev = Truncated(Normal(prev,σ_prev), 0, 1)
        prevs = rand(tn_prev, bs)
        tn_prevR = Truncated(Normal(prevR,σ_prevR), 0, 1)
        prevRs = rand(tn_prevR, bs)
        tn_mdri = Truncated(Normal(mdri,σ_mdri), 0, Inf)
        mdris = rand(tn_mdri, bs)
        tn_frr = Truncated(Normal(frr,σ_frr), 0, 1)
        frrs = rand(tn_frr, bs)
        # do stuff
        bs_incidence = kassanjee.(prevs, prevRs, mdris, frrs, T) * per
        sd = Statistics.std(bs_incidence)
        ci = Statistics.quantile(bs_incidence, [α/2, 1-α/2])
        return pe, sd, ci
    end

end

# Method for multiple surveys using one test
function incprops(prev::AbstractVector{Float64},
    prevR::AbstractVector{Float64},
    mdri::Float64, # in days
    frr::Float64; # variance-covariance matrix for prev and prevR
    σ_mdri::Float64 = 0.0,
    σ_frr::Float64 = 0.0,
    cov::Array{Float64,2} = Matrix{Float64}(I, size(prev)[1], size(prev)[1]),
    T = 730.5, # in same units as MDRI
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

# Method for multiple surveys using multiple tests
function incprops(prev::AbstractVector{Float64},
    prevR::AbstractVector{Float64},
    mdri::AbstractVector{Float64}, # in days
    frr::AbstractVector{Float64};
    σ_mdri::AbstractVector{Float64} = repeat([0.0],size(mdri)[1]),
    σ_frr::AbstractVector{Float64} = repeat([0.0],size(frr)[1]),
    cov::Array{Float64,2} = Matrix{Float64}(I, size(prev)[1], size(prev)[1]),
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

    pe, sd, ci = incprops(prev, prevR, mdri, frr, σ_prev = σ_prev,
                            σ_prevR = σ_prevR, σ_mdri = σ_mdri, σ_frr = σ_frr,
                            cov = 0.0, T = T, timeconversion = timeconversion,
                            bs = bs, α = α, per = per)
    return pe, sd, ci
end


end
