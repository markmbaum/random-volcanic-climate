using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using Statistics
using PyPlot

pygui(true)

##

V = PowerLawDistribution(7.5, 100, -1.8)
t, C, fCO2 = integrate(V)
plot(t, C .- C[1])

##

fCO2 = integrations(10_000)
println("mean = $(mean(fCO2))")
println("variance = $(var(fCO2))")
hist(𝒻T.(4.5, fCO2), bins=50, density=true);
