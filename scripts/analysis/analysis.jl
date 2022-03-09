using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using UnPack
using PyPlot

pygui(true)

##

ens = wload(datadir("sims", "ensemble.jld2"))
@unpack t, τ, σ, res = ens

##

