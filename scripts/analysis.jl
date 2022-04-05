using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using DataFrames
using Statistics
using PyPlot

pygui(true)

## load results

t, dfs = datadir("sims", "ensemble.jld2") |> loadensemble |> frameensemble
uτ = sort(unique(dfs.T.τ))
uσ = sort(unique(dfs.T.σ))

##

g = combine(groupby(dfs.T, [:τ, :σ]), :fmax => median => :fmax)
logf = log10.(reshape(g[:,:fmax], length(uτ), length(uσ)))

fig, ax = plt.subplots(1, 1, figsize=(4,3.5))
ax[:contour](
    log10.(uσ), 
    log10.(uτ), 
    logf,
    levels=[log10(3e5)],
    color="k"
)
r = ax[:pcolormesh](
    log10.(uσ), 
    log10.(uτ), 
    logf, 
    cmap="Reds",
    vmin=minimum(logf),
    vmax=6,
    shading="gouraud"
)
cb = plt.colorbar(r, ax=ax)
cb[:set_label]("log₁₀(median fCO2 peak)")
ax[:set_xlabel]("log₁₀(σ)")
ax[:set_ylabel]("log₁₀(τ)")
fig[:tight_layout]()
fig[:savefig](plotsdir("")
)

##

g = combine(
    groupby(dfs.T, [:τ, :σ]),
    :fmax => median => :fmax,
    :tsnow => mediantsnow => :tsnow
)
g = g[g.fmax .< 3e5,:]

fig, ax = plt.subplots(1,1)
cmap = plt.get_cmap("cool")
logτ = log10.(g.τ)
for h ∈ groupby(g, :τ)
    τ = h.τ[1]
    c = (log10(τ) - minimum(logτ))/(maximum(logτ) - minimum(logτ))
    ax[:semilogy](
        gya.(h.tsnow),
        h.σ,
        color=cmap(c),
        linewidth=1.5
    )
end
ax[:invert_xaxis]()
ax[:set_ylabel]("Volcanic Variance σ")
ax[:set_xlabel]("Median Time of First Snowball [Gya]")

##

figure()
cmap = plt.get_cmap("cool")
L = length(g)
for (i,k) in enumerate(g)
    semilogy(
        k.value,
        k.σ,
        label="τ=$(k.τ[1])",
        color=cmap((i-1)/(L-1)),
        linewidth=1.5
    )
end
ylabel("Volcanic Variance σ")