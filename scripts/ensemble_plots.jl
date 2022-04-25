using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using DataFrames
using Statistics
using PyPlot

pygui(true)

## load results

df = vcat(
    framevariable(:T, datadir("sims", "ensemble_wide.jld2") |> loadensemble),
    framevariable(:T, datadir("sims", "ensemble_deep.jld2") |> loadensemble)
)
uτ = sort(unique(df.τ))
uσ = sort(unique(df.σ))

##

g = combine(groupby(df, [:τ, :σ]), :fmax => median => :fmax)

fig, ax = plt.subplots(1,1)
r = ax[:contourf](
    uσ,
    uτ,
    reshape(g.fmax, length(uτ), length(uσ)),
    cmap="Oranges"
)
ax.set_yscale("log")
ax.set_xscale("log")
cb = colorbar(r, ax=ax)
cb.set_label(
    "Median CO₂ Concentration Peak [ppm]",
    rotation=270,
    va="bottom"
)
ax[:set_xlabel]("Outgassing Rate Variance (σ)")
ax[:set_ylabel]("Outgassing Rate Relaxation Time Scale, τ [yr]")
fig.tight_layout()
fig.savefig(plotsdir("fCO2_peak_map"), dpi=500)

##

ntime = size(df, 2) - findfirst(names(df) .== "1") + 1
g = combine(
    groupby(
        df,
        [:τ, :σ]
    ),
    :fmax => median => :fmax,
    [:Tmax,:Tmin] => ((a,b) -> median(a - b)) => :Trange,
    Symbol.(1:ntime) .=> (x -> quantile(x, 0.01))
)
filter!(r -> (r.fmax < 3.5e5) & (r.Trange > 2), g)

τs = [3e6, 3e7, 3e8]
fig, axs = plt.subplots(1, length(τs), figsize=(7,3.5), sharey=true)
cmap = plt.get_cmap("cool")
logσ = log10.(df.σ)
for (ax, τ) ∈ zip(axs, τs)
    sl = g[g.τ .≈ τ, :]
    for j ∈ 1:size(sl,1)
        σ = sl[j,:σ]
        c = (log10(σ) - minimum(logσ))/(maximum(logσ) - minimum(logσ))
        ax.plot(
            gya.(LinRange(2.5, 4.5, ntime)),
            values(sl[j,size(sl,2)-ntime+1:end]),
            color=cmap(c),
            linewidth=1.5,
            zorder=2
        )
    end
    ax.invert_xaxis()
    ax.set_title("τ = $(Int(τ / 1_000_000)) Myr")
    ax.set_ylim(275, 290)
end
foreach(axs) do ax
    ax.set_xlim(2, 0)
    ax.plot([2,0], [280,280], linewidth=2.5, color="k", alpha=0.5, zorder=1)
end
axs[end].annotate(raw"$T_{snow}$", (0,280.1), va="bottom", ha="right")
axs[1].set_ylabel("1ˢᵗ Temperature Percentile [K]")
fig.supxlabel("Time [Gya]")
fig.tight_layout()
cb = colorbar(
    matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.LogNorm(
            vmin=logσ |> minimum |> exp10,
            vmax=logσ |> maximum |> exp10
        ),
        cmap=cmap
    ),
    ax=axs
)
cb.set_label("Outgassing Variance (σ)")
cb.set_ticks(-5:-2 .|> exp10)
fig.savefig(plotsdir("T_quantiles"), dpi=500)
