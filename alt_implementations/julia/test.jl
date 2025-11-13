using Pkg
using Revise
Pkg.activate("./Inctools")
#push!(LOAD_PATH, "/Users/c00192/dev/inctools/julia/") # Needs an absolute path
using Inctools

prevalence(1000, 5000)
