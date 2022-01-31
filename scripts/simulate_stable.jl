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
𝒻W(::Any, t)::Float64 = -dχ(t)

#must be sure that the stepping function has no extra factors in it!
for i in 1:10
    t, C = simulate(V, 𝒻W)
    plot(t, 𝒻T.(𝒻fCO2.(C), t), "C0", alpha=0.5)
end

xlabel("Time [Gyr]")
ylabel("Temperature [K]")
