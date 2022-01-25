using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

V = PowerLawDistribution(7.5, 100, -1.8)
t, C, fCO2 = integrate(V)
plot(t, fCO2)

##

fCO2 = integrations()
hist(fCO2 .- 288);