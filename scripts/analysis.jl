using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using DataFrames
using Statistics
using PyPlot

pygui(true)

## load results

t, dfs = datadir("sims", "bigensemble.jld2") |> loadensemble |> frameensemble
uτ = sort(unique(dfs.T.τ))
uσ = sort(unique(dfs.T.σ))

##

g = groupby(
    combine(
        groupby(dfs.T, [:τ, :σ]),
        :tsnow => mediantsnow => :mediantsnow
    ),
    :τ
)

figure()
for k in g[1:4:end]
    semilogy(k.mediantsnow, k.σ, label="τ=$(k.τ[1])")
end
xlabel("Median Time of Snowball [Gya]")
ylabel("Volcanic Variance σ")
legend()