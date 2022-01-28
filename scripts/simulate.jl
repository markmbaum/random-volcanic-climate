using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

#noise term
V = PowerLaw(7.5e12, 100e12, -1.8)
#equilibrium weathering component
dχ = dΧ()

for i in 1:10
    t, C = simulate(V, (t,C) -> -dχ(t))
    plot(t, 𝒻T.(𝒻fCO2.(C), t), "C0", alpha=0.5)
end

xlabel("Time [Gyr]")
ylabel("Temperature [K]")

