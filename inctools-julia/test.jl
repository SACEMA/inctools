push!(LOAD_PATH, "/Users/eduard/dev/inctools/inctools-julia/")
#include("/Users/eduard/dev/inctools/inctools-julia/Incidence.jl")
using Incidence
using RCall

R"library(inctools)"

# DM CI
###*** DM produces too low sigma and too narrow CI!! ***###
incprops(0.2, 0.001897366596101028, 0.007, 0.0010611606852875771, 170.0, 0.1 * 170.0, 0.005, 0.5 * 0.005,
        T = 730.5, bs = 0,
        covar = 0.0,
        α = 0.05, per = 1)
R"""
incprops(PrevH = 0.2,
    RSE_PrevH = 0.001897366596101028/0.2,
    PrevR = 0.007,
    RSE_PrevR = 0.0010611606852875771/0.007,
    Boot = FALSE,
    alpha = 0.05,
    MDRI = 170,
    RSE_MDRI = 0.1,
    FRR = 0.005,
    RSE_FRR = 0.5,
    BigT = 730.5)
"""

# Bootstrapping CI
incprops(0.2, 0.001897366596101028, 0.007, 0.0010611606852875771, 170.0, 0.1 * 170.0, 0.005, 0.5 * 0.005,
        T = 730.5, bs = 100000,
        covar = 0.0,
        α = 0.05, per = 1)
R"""
incprops(PrevH = 0.2,
    RSE_PrevH = 0.001897366596101028/0.2,
    PrevR = 0.007,
    RSE_PrevR = 0.0010611606852875771/0.007,
    Boot = TRUE,
    BS_Count = 100000,
    alpha = 0.05,
    MDRI = 170,
    RSE_MDRI = 0.1,
    FRR = 0.005,
    RSE_FRR = 0.5,
    BigT = 730.5)
"""





inccounts(100000, 20000, 20000, 140, 170.0, 0.005,
    σ_mdri = 0.1 * 170.0, σ_frr = 0.5 * 0.005,
    de_npos = 1.5,
    de_nR = 1.8,
    T = 730.5, bs = 0,
    α = 0.05, per = 100)

p



p, σ = prevalence(10,100)

incprops(0.2, 0.0, 0.0068, 0.0, 170.0, 0.0, 0.005, 0.0, T = 730.5, α = 0.05)
incprops(0.3, 0.01, 170.0, 0.005, T = 730.5, α = 0.1)
incprops(0.25, 0.008, 162.0, 0.004, T = 730.5, α = 0.1)



# Try with vector of prevalences - one test
incprops([0.2, 0.3], [0.0068, 0.01], 170.0, 0.005, T = 730.5, α = 0.1)

# Try with vector of prevalences and different tests
incprops([0.2, 0.3, 0.25], [0.0068, 0.01, 0.008], [170.0,170.0,162.0], [0.005, 0.005, 0.004], T = [730.5, 730.5, 730.5], α = 0.1)

inccounts(100000, 20000, 20000, 136, 170.0, 0.005)



import Distributions
function tmvnorm(n::Int64,
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

prev = 0.2
σ_prev = 0.001897366596101028
prevR = 0.007
σ_prevR = 0.0010611606852875771
mdri = 170.0
σ_mdri = 0.1 * 170.0
frr = 0.005
σ_frr = 0.5 * 0.005
covar = 0.0
µ = [prev, prevR, mdri, frr]
Σ = [σ_prev^2 covar 0 0 ; covar σ_prevR^2 0 0 ; 0 0 σ_mdri^2 0 ; 0 0 0 σ_frr^2]

@time r = tmvnorm(1000000, µ,Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, Inf, 1.0])

vec = rand(1000000)

cov([r[:,1] vec])[1,2]

mat
r


cov([])

cov(r[:,1:2])
cor(r[:,1:2])

r = r[(r[:,1] .>= 0) .&
         (r[:,2] .>= 0) .&
         (r[:,3] .>= 0) .&
         (r[:,4] .>= 0) .&
         (r[:,1] .<= 1) .&
         (r[:,2] .<= 1) .&
         (r[:,3] .<= Inf) .&
         (r[:,4] .<= 1),:]

mu = µ
@rput mu
Sigma = Σ
@rput Sigma

R"""
library(tmvtnorm)
system.time(rtmvnorm(n = 1000000,
    mean = mu,
    sigma = Sigma,
    lower = c(0,0,0,0),
    upper = c(1, 1, Inf, 1)))
"""
