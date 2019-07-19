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
