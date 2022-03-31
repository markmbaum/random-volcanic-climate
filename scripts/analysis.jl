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

g = filter(
    x -> x.τ[1] < 5e6,
    groupby(
        combine(
            groupby(dfs.T, [:τ, :σ]),
            :tsnow => mediantsnow => :value
        ),
        :τ
    )
)

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
xlabel("Median Time of Snowball [Gya]")
ylabel("Volcanic Variance σ")
