using DrWatson
@quickactivate "Random Volcanic Climate"
using Pkg
Pkg.instantiate()

push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using IterTools: product
using Base.Threads: nthreads

#------------------------------------------------------------------------------
## INPUT/PARAMETERS

#simulation start time [Gyr]
t₁ = 2.5
#simulation end time [Gyr]
t₂ = 4.5
#values for outgassing relaxation
τ = exp10.(LinRange(6, 8, 3))
#values for outgassing variance
σ = exp10.(LinRange(-6, -4, 3))
#weathering function
𝒻W(C,t) = 𝒻whak(C, t, β=0)
#number of simulations per parameter combination
nrealize = 10*nthreads()
#number of steps for each simulation
nstep = 1_000_000
#number of time slices to store
nstore = 51

#------------------------------------------------------------------------------
## MAIN

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
