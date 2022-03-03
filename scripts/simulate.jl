using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using Statistics: mean, median
using PyPlot

pygui(true)

##

V = PowerLaw(Vᵣ, 100Vᵣ, -1.8)
𝒻W(C,t) = 𝒻whak(C, t, β=0)
t₁ = 3
C₁ = 𝒻Cₑ(t₁) #𝒻Wₑ(𝒻W, t₁, mean(V))

t, C = simulate(V, 𝒻W, C₁=C₁, t₁=t₁, nstep=100_000)
fCO2 = 𝒻fCO2.(C)
T = 𝒻T.(fCO2, t)
m = (t .> 4)
println("mean temperature = ", mean(𝒻T.(𝒻fCO2.(C),t)))
println("median temperature = ", median(𝒻T.(𝒻fCO2.(C),t)))

gya = 4.5 .- t
χ = Χ()

fig, axs = subplots(4, 1, figsize=(6,7))

axs[1][:plot](gya, C, "C0")
axs[1][:plot](gya, χ.(t), "k", alpha=0.5, zorder=2)

axs[2][:plot](gya, log10.(fCO2), "C1")
axs[2][:plot](gya, log10.(𝒻fCO2.(χ.(t))), "k", alpha=0.5, zorder=2)

axs[3][:plot](gya, T, "C2")
axs[3][:plot](gya, 𝒻T.(𝒻fCO2.(χ.(t)),t), "k", alpha=0.5, zorder=2)


axs[4][:plot](gya, 𝒻W.(C,t), ".-", color="C3")
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

figure()
hist(𝒻W.(C[end-10000:end],t[end-10000:end]), density=true, log=true, bins=40, color="gray");
