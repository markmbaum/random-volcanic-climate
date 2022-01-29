using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

#noisy CO2 outgassing
V = PowerLaw(Vᵣ, 10Vᵣ, -1.8)
#equilibrium weathering component
dχ = dΧ()
#actual weathering function
𝒻W(t, ::Any)::Float64 = -dχ(t)

for i in 1:10
    t, C = simulate(V, 𝒻W)
    plot(t, 𝒻T.(𝒻fCO2.(C), t), "C0", alpha=0.5)
end

xlabel("Time [Gyr]")
ylabel("Temperature [K]")
