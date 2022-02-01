using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

V = PowerLaw(Vᵣ, 5000Vᵣ, -1.8)
𝒻W(C,t) = 𝒻whak(C, t, β=0)
t₁ = 2.5
C₁ = 𝒻Wₑ(𝒻W, t₁, mean(V))

fig, axs = subplots(4, 1, figsize=(6,8))
t, C = simulate(V, 𝒻W, C₁=C₁, t₁=t₁)

gya = 4.5 .- t
χ = Χ()

axs[1][:plot](gya, C, "C0")
axs[1][:plot](gya, χ.(t), "k", alpha=0.5, zorder=2)

axs[2][:plot](gya, log10.(𝒻fCO2.(C)), "C1")
axs[2][:plot](gya, log10.(𝒻fCO2.(χ.(t))), "k", alpha=0.5, zorder=2)

axs[3][:plot](gya, 𝒻T.(𝒻fCO2.(C),t), "C2")
axs[3][:plot](gya, 𝒻T.(𝒻fCO2.(χ.(t)),t), "k", alpha=0.5, zorder=2)


axs[4][:plot](gya, 𝒻W.(C,t), "C3")
axs[4][:plot](gya, fill(mean(V), length(t)), "k", alpha=0.5, zorder=2)

axs[1][:set_title]("Ocean-Atm Carbon [Tmole]")
axs[2][:set_title]("log₁₀ fCO2 [ppmv]")
axs[3][:set_title]("Temperature [K]")
axs[4][:set_title]("Weathering [Tmole/yr]")
axs[4][:set_xlabel]("Time [Gyr]")
for ax ∈ axs
    ax[:invert_xaxis]()
end
fig[:tight_layout]()
