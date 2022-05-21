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
        nstep=50_000
    )
    return t, V
end

##

τ = [3e6, 3e7, 3e8]
σ = [5e-4, 3e-4, 1e-4]
nstack = 5
fig, axs = plt.subplots(1, length(τ), figsize=(7,3.5), sharey=true, constrained_layout=true)
for (i,ax) in enumerate(axs)
    for j in 1:nstack
        t, V = outgassing(τ[i], σ[i])
        ax.plot(𝒻gya.(t), V .+ 6*(j-1), color="C3", linewidth=0.75, alpha=0.9)
        ax.set_title("τ = $(Int(τ[i] / 1_000_000)) Myr\nσ = $(σ[i])")
    end
end
axs[1].set_yticks([])
axs[1].set_ylabel("CO₂ Outgassing Rate [-]")
fig.supxlabel("Time [Gya]")
fig.savefig(plotsdir("outgassing_examples.png"), dpi=500)