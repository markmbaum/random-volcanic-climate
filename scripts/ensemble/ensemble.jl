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

function prettylogrange(o₁::Int, o₂::Int)
    #values within each order of magnitude
    x = [1, 1.5, 2, 2.5, 3, 4, 5, 6, 8]
    #fill in orders
    v = vcat(map(i -> x * (10. ^ i), o₁:o₂)...)
    append!(v, 10. ^ (o₂ + 1))
    return v
end

##-----------------------------------------------------------------------------
# ENSEMBLE PARAMETERS

#simulation start time [Gyr]
t₁ = 2.5
#simulation end time [Gyr]
t₂ = 4.5
#values for outgassing relaxation
τ = [1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9]
#values for outgassing variance
σ = prettylogrange(-5, -3)
#weathering function
𝒻W(C,t) = 𝒻whak(C, t, β=0)
#number of simulations per parameter combination
nrealize = 2000*nthreads()
#number of steps for each simulation
nstep = 1_000_000
#number of time slices to store
nstore = 17

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
