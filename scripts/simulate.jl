using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

V = PowerLaw(Vᵣ/50, 10Vᵣ, -1.8)
𝒻W(C,t) = 𝒻whak(C, t)
t₁ = 2.5
C₁ = 𝒻Cₑ(t₁)

fig, axs = subplots(3, 1, figsize=(6,7))
t, C = simulate(V, 𝒻W, C₁=C₁, t₁=t₁)
gya = 4.5 .- t
axs[1][:plot](gya, log10.(𝒻fCO2.(C)), "C0")
axs[2][:plot](gya, 𝒻T.(𝒻fCO2.(C),t), "C1")
axs[3][:plot](gya, 𝒻W.(C,t), "C2")
axs[3][:plot](
    [minimum(gya), maximum(gya)],
    [mean(V), mean(V)],
    "r--",
    alpha=0.8,
    zorder=-1,
    label="Outgassing"
)
axs[1][:set_ylabel]("log₁₀ fCO2 [ppmv]")
axs[2][:set_ylabel]("Temperature [K]")
axs[3][:set_ylabel]("Weathering [teramole/yr]")
axs[3][:set_xlabel]("Time [Gyr]")
for ax ∈ axs
    ax[:invert_xaxis]()
end
fig[:tight_layout]()
