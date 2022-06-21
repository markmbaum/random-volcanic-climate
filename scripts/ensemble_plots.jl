using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using DataFrames
using Statistics
using PyPlot

pygui(true)

## load results

enswide = datadir("sims", "ensemble_wide.jld2") |> loadensemble
ensdeep = datadir("sims", "ensemble_deep.jld2") |> loadensemble
T = vcat(
    framevariable(:T, enswide),
    framevariable(:T, ensdeep)
)
Tₑ = vcat(
    framevariable(wload(datadir("exp_pro", "Te_wide.jld2"))["Tₑ"], enswide),
    framevariable(wload(datadir("exp_pro", "Te_deep.jld2"))["Tₑ"], ensdeep)
)
uτ = T.τ |> unique |> sort
uσ = T.σ |> unique |> sort
@assert all(enswide[:t] .== ensdeep[:t])
time = enswide[:t]
tc = timecols(T)

##

g = combine(groupby(T, [:τ, :σ]), :fmax => median => :fmax)

fig, ax = plt.subplots(1, 1, figsize=(5,5), constrained_layout=true)
r = ax[:pcolormesh](
    uσ,
    uτ/1e6,
    reshape(g.fmax, length(uτ), length(uσ)),
    cmap="Oranges",
    shading="gouraud",
    vmax=1e6,
    vmin=0
)
ax.set_yscale("log")
ax.set_xscale("log")
cb = colorbar(r, ax=ax)
ax.set_title("Median CO₂ Concentration Peak [ppm]")
ax.set_xlabel("Outgassing Rate SD, σ [Tmole/yr]")
ax.set_ylabel("Outgassing Rate Relaxation Time Scale, τ [Myr]")
fig.savefig(plotsdir("fCO2_peak_map"), dpi=500)

##

g = combine(
    groupby(
        T,
        [:τ, :σ]
    ),
    :fmax => median => :fmax,
    [:Tmax,:Tmin] => ((a,b) -> median(a - b)) => :Trange,
    tc .=> (x -> quantile(x, 0.01)) .=> tc
)
filter!(r -> (r.fmax < 4.2e5) & (r.Trange > 1), g)

τs = [3e6, 1e7, 3e7, 1e8, 3e8]
fig, axs = plt.subplots(1, length(τs), figsize=(7,3.5), sharey=true, constrained_layout=true)
cmap = plt.get_cmap("cool")
logσ = log10.(T.σ)
𝒻color(σ) = cmap((log10(σ) - minimum(logσ))/(maximum(logσ) - minimum(logσ)))
for (ax, τ) ∈ zip(axs, τs)
    sl = g[g.τ .≈ τ, :]
    for j ∈ 1:size(sl,1)
        σ = sl[j,:σ]
        ax.plot(
            𝒻gya.(time),
            values(sl[j,size(sl,2)-length(time)+1:end]),
            color=𝒻color(σ),
            linewidth=1.5,
            zorder=2
        )
    end
    ax.invert_xaxis()
    ax.set_title("τ = $(Int(τ / 1_000_000)) Myr", fontsize=10)
    ax.set_ylim(275, 290)
end
foreach(axs) do ax
    ax.set_xlim(2, 0)
    ax.set_xticks([2,0])
    ax.plot([2,0], [280,280], linewidth=2.5, color="k", alpha=0.5, zorder=1)
end
axs[end].annotate(
    raw"$T_{snow}$",
    (0,279.8),
    va="top",
    ha="right"
)
axs[1].set_ylabel("1ˢᵗ Temperature Percentile [K]")
fig.supxlabel("Time [Gya]")
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
cb.set_label(
    "Outgassing Rate std, σ [Tmole/yr]",
    rotation=270,
    va="bottom"
)
cb.set_ticks(-5:-2 .|> exp10)
fig.savefig(plotsdir("T_quantiles"), dpi=500)

##
