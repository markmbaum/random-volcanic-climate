using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using DataFrames
using PyPlot

pygui(true)

##

t, dfC, dfV = loadensemble(
    datadir(
        "sims",
        "ensemble.jld2"
    )
);

##
