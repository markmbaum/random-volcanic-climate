using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

#noisy CO2 outgassing
V = PowerLaw(Vᵣ, 10Vᵣ, -1.8)

fig, axs = subplots(3, 1, figsize=(6,7))
for i in 1:4
    t, C = simulate(V, 𝒻mac)
    axs[1][:plot](4.5 .- t, log10.(𝒻fCO2.(C)), "C0", alpha=0.5)
    axs[2][:plot](4.5 .- t, 𝒻T.(𝒻fCO2.(C),t), "C1", alpha=0.5)
    axs[3][:plot](4.5 .- t, 𝒻mac.(t, C), "C2", alpha=0.5)
    axs[3][:plot](4.5 .- [t[1], t[end]], [mean(V), mean(V)], "k", alpha=0.3)
end
axs[1][:set_ylabel]("log₁₀ fCO2 [ppmv]")
axs[2][:set_ylabel]("Temperature [K]")
axs[3][:set_ylabel]("Weathering [teramole/yr]")
axs[3][:set_xlabel]("Time [Gyr]")
for ax ∈ axs
    ax[:invert_xaxis]()
end
fig[:tight_layout]()
