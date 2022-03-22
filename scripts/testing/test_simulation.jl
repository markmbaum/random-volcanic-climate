using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using Statistics: mean, median
using PyPlot

pygui(true)

##

𝒻W(C,t) = 𝒻whak(C, t, β=0)
t₁ = 2.5

t, C, V = simulate(
    initparams(
        μ=Vᵣ,
        τ=1e7,
        σ=1e-4,
        Vₘ=0
    ),
    t₁=t₁,
    C₁=𝒻Cₑ(t₁),
    𝒻W=𝒻W,
    nstep=100_000
)
println(minimum(V))

fCO2 = 𝒻fCO2.(C)
T = 𝒻T.(fCO2, t)
#m = (t .> 4)
#println("mean temperature = ", mean(𝒻T.(𝒻fCO2.(C),t)))
#println("median temperature = ", median(𝒻T.(𝒻fCO2.(C),t)))

gya = 4.5 .- t
χ = Χ()

fig, axs = subplots(5, 1, figsize=(6,7))

axs[1][:plot](gya, V, "C3")
axs[1][:plot]([gya[1], gya[end]], [Vᵣ, Vᵣ], "k", alpha=0.5, zorder=2)

axs[2][:plot](gya, C, "C0")
axs[2][:plot](gya, χ.(t), "k", alpha=0.5, zorder=2)

axs[3][:plot](gya, log10.(fCO2), "C1")
axs[3][:plot](gya, log10.(𝒻fCO2.(χ.(t))), "k", alpha=0.5, zorder=2)

axs[4][:plot](gya, T, "C2")
axs[4][:plot]([gya[1], gya[end]], [Tᵣ, Tᵣ], "k", alpha=0.5, zorder=2)

axs[5][:plot](gya, 𝒻W.(C,t), "-", color="C4")
axs[5][:plot](gya, fill(mean(V), length(t)), "k", alpha=0.5, zorder=2)

axs[1][:set_title]("Outgassing Rate [Tmole/yr]")
axs[2][:set_title]("Ocean-Atm Carbon [Tmole]")
axs[3][:set_title]("log₁₀ fCO2 [ppmv]")
axs[4][:set_title]("Temperature [K]")
axs[5][:set_title]("Weathering [Tmole/yr]")
axs[5][:set_xlabel]("Time [Gyr]")
for ax ∈ axs
    ax[:invert_xaxis]()
end
fig[:tight_layout]()

##

figure()
hist(T, density=true, log=true, bins=40, color="gray");
