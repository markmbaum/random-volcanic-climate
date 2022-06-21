using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

t₁ = 2.5
𝒻W(C,t) = 𝒻whak(C, t, β=0)

t, C, V = simulate(
    initparams(
        μ=Vᵣ,
        τ=30e6,
        σ=2e-4,
        Vₘ=0
    ),
    t₁=t₁,
    C₁=𝒻Cₑ(t₁),
    𝒻W=𝒻W,
    nstep=1_000_000
)

fCO2 = 𝒻fCO2.(C)
T = 𝒻T.(fCO2, t)

gya = 𝒻gya.(t)
χ = Χ()

fig, axs = plt.subplots(5, 1, figsize=(5.5,7), constrained_layout=true)

axs[1].plot(gya, V, "C3")
axs[1].plot([gya[1], gya[end]], [Vᵣ, Vᵣ], "k", alpha=0.5, zorder=2)

axs[2].plot(gya, log10.(C), "C0")
axs[2].plot(gya, log10.(χ.(t)), "k", alpha=0.5, zorder=2)

axs[3].plot(gya, log10.(fCO2), "C1")
axs[3].plot(gya, log10.(𝒻fCO2.(χ.(t))), "k", alpha=0.5, zorder=2)

axs[4].plot(gya, T, "C2")
axs[4].plot([gya[1], gya[end]], [Tᵣ, Tᵣ], "k", alpha=0.5, zorder=2)

axs[5].plot(gya, 𝒻W.(C,t), "-", color="C4")
axs[5].plot(gya, fill(Vᵣ, length(t)), "k", alpha=0.5, zorder=2)

axs[1][:set_title]("Outgassing Rate [Tmole/yr]")
axs[2][:set_title]("log₁₀ Ocean-Atmosphere Carbon [Tmole]")
axs[3][:set_title]("log₁₀ fCO2 [ppmv]")
axs[4][:set_title]("Temperature [K]")
axs[5][:set_title]("Weathering [Tmole/yr]")
axs[5][:set_xlabel]("Time [Gyr]")
for (ax,c) ∈ zip(axs, "abcde")
    ax.annotate(
        string(c),
        (0.014,1.05),
        xycoords="axes fraction",
        va="top",
        fontsize=12,
        fontweight="bold",
        backgroundcolor="white"
    )
    ax.invert_xaxis()
    ax.set_xlim(2, 0)
end
fig.savefig(plotsdir("example_simulation"), dpi=500)
plt.close(fig)

##

plt.figure()
hist((@. fCO2 - 𝒻fCO2(χ(t))), density=true, log=true, bins=40, color="gray");
