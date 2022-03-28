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
gya = 4.5 .- t

## quantiles over time for each τ,σ combination

uτ = sort(unique(dfs.T.τ))
uσ = sort(unique(dfs.T.σ))
Q = [0.001, 0.01, 0.1, 0.25]

fig, axs = plt.subplots(length(uτ), length(uσ), sharey=true)
for (i,τ) ∈ enumerate(uτ), (j,σ) ∈ enumerate(uσ)
    sl = filter(r -> (r.τ ≈ τ) & (r.σ ≈ σ), dfs.T) 
    for q ∈ Q
        axs[i,j].plot(
            gya,
            [quantile(sl[:,k+3], q) for k ∈ 1:length(t)]
        )
        axs[i,j].plot(
            gya,
            [quantile(sl[:,k+3], 1-q) for k ∈ 1:length(t)]
        )
    end
    axs[i,j].invert_xaxis()
    if i < length(uτ)
        axs[i,j].set_xticklabels([])
    end
    if j == 1
        axs[i,j].set_ylabel("τ = $τ")
    end
    if i == length(uτ)
        axs[i,j].set_xlabel("σ = $σ")
    end
    axs[i,j].set_title(round(Float64(τ*σ^2), sigdigits=4))
end
fig.suptitle("Temperature Quantiles over Time")
fig.tight_layout()

##

dgs = [x[:,["τ","σ","1",string(length(t))]] for x in groupby(dfs.T, [:τ, :σ])]

fig, axs = plt.subplots(1, length(Q))
for (i,q) ∈ enumerate(Q)
    ax = axs[i]    
    x = [log10(dg.σ[1]/dg.τ[1]) for dg ∈ dgs]
    y = [std(dg[:,3]) - std(dg[:,4]) for dg in dgs]
    ax.scatter(x, y)
end
fig.tight_layout()
