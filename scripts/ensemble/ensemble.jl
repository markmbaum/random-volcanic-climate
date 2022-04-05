using DrWatson
@quickactivate "Random Volcanic Climate"
using Pkg
Pkg.instantiate() 

push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using IterTools: product
using Base.Threads: nthreads

##-----------------------------------------------------------------------------
# FUNCTIONS

function prettylogrange(o₁, o₂)
    #values within each order of magnitude
    x = [1, 1.5, 2, 2.5, 3, 4, 5, 6, 8]
    #fill in orders
    vcat(map(i -> x * (10. ^ i), o₁:o₂)...)
end

##-----------------------------------------------------------------------------
# ENSEMBLE PARAMETERS

#simulation start time [Gyr]
t₁ = 2.5
#simulation end time [Gyr]
t₂ = 4.5
#values for outgassing relaxation
τ = prettylogrange(5, 8)
#values for outgassing variance
σ = exp10.(LinRange(-6, -3, 100)) #prettylogrange(-6, -3)
#weathering function
𝒻W(C,t) = 𝒻whak(C, t, β=0)
#number of simulations per parameter combination
nrealize = 50*nthreads()
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
