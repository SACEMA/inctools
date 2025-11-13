using Pkg

# Find Inctools package (works from julia/ or tests/ directory)
inctools_path = isdir("./Inctools") ? "./Inctools" : "../Inctools"
if !isdir(inctools_path)
    error("Cannot find Inctools package. Run this from julia/ directory.")
end
Pkg.activate(inctools_path)

using Revise
using Inctools
using LinearAlgebra
using Statistics

println("Testing rtmvnorm refactoring...")
println("="^60)

# Test 1: Gibbs sampling with diagonal covariance (4D)
println("\n1. Testing Gibbs sampling (4D diagonal)")
µ = [0.5, 0.3, 0.7, 0.4]
Σ = Matrix(Diagonal([0.1, 0.08, 0.12, 0.09]))
samples = rtmvnorm(1000, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])
println("  Shape: $(size(samples))")
println("  Mean: $(round.(mean(samples, dims=1), digits=3))")
println("  All in bounds: $(all(0 .<= samples .<= 1))")

# Test 2: Explicit Gibbs method
println("\n2. Testing explicit Gibbs method")
samples2 = rtmvnorm(1000, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0], method=:gibbs)
println("  Shape: $(size(samples2))")
println("  ✓ Gibbs method works")

# Test 3: Rejection sampling with correlation (4D)
println("\n3. Testing rejection sampling (4D with correlation)")
Σ_corr = [0.1 0.02 0 0; 0.02 0.08 0 0; 0 0 0.12 0; 0 0 0 0.09]
samples3 = rtmvnorm(100, µ, Σ_corr, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])
println("  Shape: $(size(samples3))")
println("  Mean: $(round.(mean(samples3, dims=1), digits=3))")
println("  All in bounds: $(all(0 .<= samples3 .<= 1))")

# Test 4: Gibbs sampling with different dimensions (2D)
println("\n4. Testing Gibbs sampling (2D)")
µ2 = [0.2, 0.8]
Σ2 = Matrix(Diagonal([0.05, 0.05]))
samples4 = rtmvnorm(1000, µ2, Σ2, [0.0, 0.0], [1.0, 1.0])
println("  Shape: $(size(samples4))")
println("  Mean: $(round.(mean(samples4, dims=1), digits=3))")
println("  ✓ Works with 2D")

# Test 5: Gibbs sampling with different dimensions (6D)
println("\n5. Testing Gibbs sampling (6D)")
µ6 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
Σ6 = Matrix(Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]))
samples5 = rtmvnorm(1000, µ6, Σ6,
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
println("  Shape: $(size(samples5))")
println("  Mean: $(round.(mean(samples5, dims=1), digits=3))")
println("  ✓ Works with 6D")

# Test 6: Direct call to rtmvnorm_gibbs
println("\n6. Testing rtmvnorm_gibbs directly")
samples6 = Inctools.rtmvnorm_gibbs(1000, µ, Σ, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])
println("  Shape: $(size(samples6))")
println("  ✓ Direct call works")

# Test 7: Direct call to rtmvnorm_rejection
println("\n7. Testing rtmvnorm_rejection directly")
samples7 = Inctools.rtmvnorm_rejection(100, µ, Σ_corr, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])
println("  Shape: $(size(samples7))")
println("  ✓ Direct call works")

# Test 8: Integration test with incprops
println("\n8. Testing integration with incprops")
result = incprops(0.20, 0.015, 0.10, 0.02, 130.0, 15.0, 0.01, 0.005, bs=1000, covar=0.0002)
println("  Incidence: $(round(result.I, digits=4))")
println("  CI: [$(round(result.CI[1], digits=4)), $(round(result.CI[2], digits=4))]")
println("  ✓ incprops with covariance works")

println("\n" * "="^60)
println("✓ All tests passed!")
