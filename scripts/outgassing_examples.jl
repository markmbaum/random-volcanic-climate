using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

function outgassing(τ, σ)
    t, C, V = simulate(
        initparams(
            μ=Vᵣ,
            τ=τ,
            σ=σ,
            Vₘ=0
        ),
        t₁=2.5,
        C₁=𝒻Cₑ(2.5),
        𝒻W=(C,t)->𝒻whak(C, t, β=0),
        nstep=10_000
    )
    return t, V
end

##

fig, axs = plt.subplots(1, 3, figsize=(7,3.5), sharey=true)
τ = [3e6, 3e7, 3e8]
σ = [5e-4, 2.5e-4, 1e-4]
for (i,ax) in enumerate(axs)
    t, V = outgassing(τ[i], σ[i])
    ax.plot(𝒻gya.(t), V, color="C3", linewidth=0.75, alpha=0.7)
    ax.set_title("τ = $(Int(τ[i] / 1_000_000)) Myr\nσ = $(σ[i])")
end
axs[1].set_ylabel("CO₂ Outgassing Rate [Tmole/yr]")
fig.supxlabel("Time [Gya]")
fig.tight_layout()
#fig.savefig(plotsdir("outgassing_examples.png"), dpi=500)