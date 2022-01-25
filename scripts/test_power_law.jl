using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

## check normalization

for μ ∈ 1:0.5:5, x₂ ∈ 6:0.5:10, γ ∈ -1.1:-0.2:-3.1
    D = PowerLawDistribution(μ, x₂, γ)
    q = ∫(D)
    if !(q ≈ 1)
        error("not normalized with q=$q μ=$μ x₂=$x₂ γ=$γ")
    end
end
println("normalization is working as expected")

## inspect a distribution

μ = 5
x₂ = 100
γ = -1.8

D = PowerLawDistribution(μ, x₂, γ)
hist(D(1_000_000), bins=50, density=true, log=true);
