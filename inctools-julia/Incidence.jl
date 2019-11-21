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

module Incidence

export prevalence, rtmvnorm, incprops, inccounts, incdif

import LinearAlgebra.I
import Distributions
import Statistics
import DataFrames

function prevalence(pos, n, de = 1; ci = false, α = 0.05) #, f = 1
    if n == 0
        @warn "n = 0: prevalence undefined. n set to 1."
        n = 1
    end
    if pos < 5
        @warn "Too few successes for the normal approximation to be valid"
    end
    if n - pos < 5
        @warn "Sample size twoo small for the normal apprximation to be valid"
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

# This function does not implement the appropriate methd
# !! Replace with a gibbs sampler!
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
            @warn "Truncated multivariate normal with covariance not implemented yet"
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
                    rand(Distributions.TruncatedNormal(prev, σ_prev, 0.0, 1.0), bs),
                    rand(Distributions.TruncatedNormal(prevR, σ_prevR, 0.0, 1.0), bs),
                    rand(Distributions.TruncatedNormal(mdri, σ_mdri, 0.0, Inf), bs),
                    rand(Distributions.TruncatedNormal(frr, σ_frr, 0.0, 1.0), bs)
                    )
        elseif covar > 0.0
            @warn "Truncated multivariate normal with covariance not implemented yet"
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
        dmdri = Distributions.TruncatedNormal(mdri, σ_mdri, 0, Inf)
        dfrr = Distributions.TruncatedNormal(frr, σ_frr, 0, 1)
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
                    rand(Distributions.TruncatedNormal(prev[1], σ_prev[1], 0.0, 1.0), bs),
                    rand(Distributions.TruncatedNormal(prevR[1], σ_prevR[1], 0.0, 1.0), bs),
                    rand(Distributions.TruncatedNormal(prev[2], σ_prev[2], 0.0, 1.0), bs),
                    rand(Distributions.TruncatedNormal(prevR[2], σ_prevR[2], 0.0, 1.0), bs),
                    rand(Distributions.TruncatedNormal(mdri, σ_mdri, 0.0, Inf), bs),
                    rand(Distributions.TruncatedNormal(frr, σ_frr, 0.0, 1.0), bs)
                    )
        elseif covar[1] > 0.0 || covar[2] > 0.0
            @warn "Truncated multivariate normal with covariance not implemented yet"
            µ = [prev[1], prevR[1], prev[2], prevR[2], mdri, frr]
            Σ = [σ_prev[1]^2 covar[1] 0 0 0 0; covar[1] σ_prevR[1]^2 0 0 0 0 ; 0 0 σ_prev[2]^2 covar[2] 0 0 ; 0 0 covar[2] σ_prevR[2]^2 0 0 ; 0 0 0 0 σ_mdri^2 0 ; 0 0 0 0 0 σ_frr^2]
            r = rtmvnorm(bs, µ, Σ, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0, Inf, 1.0])
        end

        bs_incidence_1 = kassanjee.(r[:,1], r[:,2], r[:,5], r[:,6], T) .* per
        bs_incidence_2 = kassanjee.(r[:,3], r[:,4], r[:,5], r[:,6], T) .* per
        if any(x-> x < 0, bs_incidence_1) || any(x-> x < 0, bs_incidence_2)
            @warn "Negative boostrapped incidence estimates set to 0.0 for difference calculations"
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
        dmdri = Distributions.TruncatedNormal(mdri, σ_mdri, 0, Inf)
        dfrr = Distributions.TruncatedNormal(frr, σ_frr, 0, 1)
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
            @warn "Negative boostrapped incidence estimates set to 0.0 for difference calculations"
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
