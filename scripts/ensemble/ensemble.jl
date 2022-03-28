using DrWatson
@quickactivate "Random Volcanic Climate"
using Pkg
Pkg.instantiate() 

push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using IterTools: product
using Base.Threads: nthreads

##-----------------------------------------------------------------------------
# INPUT/PARAMETERS

#simulation start time [Gyr]
t₁ = 2.5
#simulation end time [Gyr]
t₂ = 4.5
#values for outgassing relaxation
τ = [
    1e5,
    2e5,
    4e5,
    1e6,
    2e6,
    2.5e6,
    4e6,
    1e7,
    2e7,
    4e7,
    1e8,
    2e8,
    4e8
]
#values for outgassing variance
σ = [
    1e-6,
    2e-6,
    4e-6,
    1e-5,
    2e-5,
    4e-5,
    1e-4,
    2e-4,
    4e-4,
    8e-4,
    1e-3,
    2e-3,
    3e-3,
    4e-3,
    5e-3,
    6e-3,
    8e-3,
    1e-2
]
#weathering function
𝒻W(C,t) = 𝒻whak(C, t, β=0)
#number of simulations per parameter combination
nrealize = 500*nthreads()
#number of steps for each simulation
nstep = 1_000_000
#number of time slices to store
nstore = 11

##-----------------------------------------------------------------------------
# MAIN

#create parameter combinations
params = product(τ, σ)

##

#simulate ensemble and save directly
saveensemble(
    datadir(
        "sims",
        "ensemble.jld2"
    ),
    ensemble(
        params,
        t₁,
        t₂,
        nrealize,
        nstep,
        nstore,
        𝒻W
    )
)
