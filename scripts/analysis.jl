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

df = combine(
    groupby(
        dfs.T,
        [:τ, :σ]
    ),
    :fmax => median => :fmax,
    :Tmax => median => :Tmax,
    :Tmin => median => :Tmin
)
f = reshape(df[:,:fmax], length(uτ), length(uσ))
Tmax = reshape(df[:,:Tmax], length(uτ), length(uσ))
Tmin = reshape(df[:,:Tmin], length(uτ), length(uσ))

fig, axs = plt.subplots(3, 1, figsize=(3.25,9))

r = axs[1][:pcolormesh](
    log10.(uσ), 
    log10.(uτ), 
    f,
    vmin=minimum(f),
    vmax=1e6,
    cmap="Oranges",
    shading="gouraud"
)
cb = plt.colorbar(r, ax=axs[1])
cb[:set_label]("median fCO2 peak [ppm]")
cb[:set_ticks](250e3:125e3:1e6)
axs[1][:set_ylabel]("log₁₀(τ)")

r = axs[2][:pcolormesh](
    log10.(uσ), 
    log10.(uτ), 
    Tmin,
    #vmin=minimum(Tmin),
    #vmax=maximum(Tmax),
    cmap="Blues_r",
    shading="gouraud"
)
axs[2][:set_ylabel]("log₁₀(τ)")
cb = plt.colorbar(r, ax=axs[2])
cb[:set_label]("median Temperature low [K]")

r = axs[3][:pcolormesh](
    log10.(uσ), 
    log10.(uτ), 
    Tmax,
    #vmin=minimum(Tmin),
    #vmax=maximum(Tmax),
    cmap="Reds",
    shading="gouraud"
)
axs[3][:set_ylabel]("log₁₀(τ)")
axs[3][:set_xlabel]("log₁₀(σ)")
cb = plt.colorbar(r, ax=axs[3])
cb[:set_label]("median temperature peak [ppm]")

fig[:tight_layout]()
#fig[:savefig](plotsdir("fCO2_peaks"), dpi=500)

##

df = combine(
    groupby(
        dfs.T,
        [:τ, :σ]
    ),
    :fmax => median => :fmax,
    :tsnow => mediantsnow => :tsnow
)
df = df[df.fmax .< 4e5,:]

fig, ax = plt.subplots(1,1)
cmap = plt.get_cmap("cool")
logτ = log10.(df.τ)
for h ∈ groupby(df, :τ)
    τ = h.τ[1]
    c = (log10(τ) - minimum(logτ))/(maximum(logτ) - minimum(logτ))
    ax[:semilogy](
        gya.(h.tsnow),
        h.σ,
        color=cmap(c),
        linewidth=1.5
    )
end
cb = colorbar(
    matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.LogNorm(
            vmin=logτ |> minimum |> exp10,
            vmax=logτ |> maximum |> exp10
        ),
        cmap=cmap
    ),
    ax=ax
)
cb.set_label("Outgassing Relaxation Time Scale (τ)")
ax[:invert_xaxis]()
ax[:set_ylabel]("Volcanic Variance (σ)")
ax[:set_xlabel]("Median Time of First Snowball [Gya]")
fig.tight_layout()
fig.savefig(plotsdir("median_tsnow"), dpi=500)

##

ntime = size(dfs.T, 2) - 6
df = combine(
    groupby(
        dfs.T, 
        [:τ, :σ]
    ),
    :fmax => median => :fmax,
    Symbol.(1:ntime) .=> (x -> quantile(x, 0.005))
)
df = df[df.fmax .< 3.5e5, Not(:fmax)]

τs = [1e6, 3e6, 1e7, 3e7, 1e8]
fig, axs = plt.subplots(1, length(τs), figsize=(7,3), sharey=true)
cmap = plt.get_cmap("cool")
logσ = log10.(df.σ)
for (ax, τ) ∈ zip(axs, τs)
    sl = df[df.τ .≈ τ, :]
    for j ∈ 1:size(sl,1)
        σ = sl[j,:σ]
        c = (log10(σ) - minimum(logσ))/(maximum(logσ) - minimum(logσ))
        ax.plot(
            gya.(LinRange(2.5, 4.5, ntime)),
            values(sl[j,3:end]),
            color=cmap(c),
            linewidth=1.5
        )
    end
    ax.invert_xaxis()
    ax.set_title("τ = $(τ / 1_000_000)\nMyr")
    ax.set_ylim(275, 290)
end
axs[1].set_ylabel("Temperature [K]")
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
fig.savefig(plotsdir("T_quantiles"), dpi=500)

##

ntime = size(dfs.T, 2) - 6
df = combine(
    groupby(
        dfs.T, 
        [:τ, :σ]
    ),
    :fmax => median => :fmax,
    Symbol.(1:ntime) .=> (x -> count(y -> y < 280, x)/length(x))
)
df = df[df.fmax .< 4e5, Not(:fmax)]

τs = [1e6, 3e6, 1e7, 3e7, 1e8]
fig, axs = plt.subplots(1, length(τs), figsize=(7,3), sharey=true)
cmap = plt.get_cmap("cool")
logσ = log10.(df.σ)
for (ax, τ) ∈ zip(axs, τs)
    sl = df[df.τ .≈ τ, :]
    for j ∈ 1:size(sl,1)
        σ = sl[j,:σ]
        c = (log10(σ) - minimum(logσ))/(maximum(logσ) - minimum(logσ))
        ax.semilogy(
            gya.(LinRange(2.5, 4.5, ntime)),
            values(sl[j,3:end]),
            color=cmap(c),
            linewidth=1.5
        )
    end
    ax.invert_xaxis()
    ax.set_title("τ = $(τ / 1_000_000)\nMyr")
end
axs[1].set_ylabel("Snowball Fraction")
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
fig.savefig(plotsdir("snowball_prob"), dpi=500)
