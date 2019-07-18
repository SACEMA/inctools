push!(LOAD_PATH, "/Users/eduard/dev/inctools/inctools-julia/")
#include("/Users/eduard/dev/inctools/inctools-julia/Incidence.jl")
using Incidence

p, σ = prevalence(10,100)

incprops(0.2, 0.0068, 170.0, 0.005, T = 730.5, α = 0.1)
incprops(0.3, 0.01, 170.0, 0.005, T = 730.5, α = 0.1)
incprops(0.25, 0.008, 162.0, 0.004, T = 730.5, α = 0.1)

# Try with vector of prevalences - one test
incprops([0.2, 0.3], [0.0068, 0.01], 170.0, 0.005, T = 730.5, α = 0.1)

# Try with vector of prevalences and different tests
incprops([0.2, 0.3, 0.25], [0.0068, 0.01, 0.008], [170.0,170.0,162.0], [0.005, 0.005, 0.004], T = [730.5, 730.5, 730.5], α = 0.1)




inccounts(100000, 20000, 20000, 136, 170.0, 0.005)


# bootstap for CI
incprops(0.2, 0.007, 170.0, 0.005,
    σ_prev = 0.001897366596101028, σ_prevR = 0.0010611606852875771,
    σ_mdri = 0.1 * 170.0, σ_frr = 0.5 * 0.005,
    T = 730.5, bs = 100000,
    α = 0.05, per = 100)

incprops(0.2, 0.007, 170.0, 0.005,
        σ_prev = 0.001897366596101028, σ_prevR = 0.0010611606852875771,
        σ_mdri = 0.1 * 170.0, σ_frr = 0.5 * 0.005,
        T = 730.5, bs = 0,
        α = 0.05, per = 100)

inccounts(100000, 20000, 20000, 140, 170.0, 0.005,
    σ_mdri = 0.1 * 170.0, σ_frr = 0.5 * 0.005,
    de_npos = 1.5,
    de_nR = 1.8,
    T = 730.5, bs = 0,
    α = 0.05, per = 100)

prevalence(2000,10000)
prevalence(14,2000)
