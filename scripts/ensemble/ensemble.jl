using DrWatson
@quickactivate "Random Volcanic Climate"

using Pkg
Pkg.instantiate() 

push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using DataFrames
using Statistics
using IterTools: product
using Base.Threads: nthreads

##-----------------------------------------------------------------------------
# CONSTANTS/INPUTS

#simulation start time [Gyr]
const t₁ = 2.5
#simulation end time [Gyr]
const t₂ = 4.5
#weathering function
const 𝒻W(C,t) = 𝒻whak(C, t, β=0)
#number of steps for each simulation
const nstep = 1_000_000
#number of time slices to store
const nstore = 17
#values for outgassing relaxation
const τ = prettylogrange(5, 8, [1, 2.15, 3, 4.6])
#values for outgassing variance
const σ = prettylogrange(-5, -3, [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8])
#number of realizations per parameter combo in initial/wide ensemble
nwide = 60*nthreads()
#number of realizations per parameter combo in second/deep ensemble
ndeep = 2000*nthreads()

##-----------------------------------------------------------------------------
# FUNCTIONS

function prettylogrange(o₁::Int, o₂::Int, x::AbstractArray)
    v = vcat(map(i -> x * (10. ^ i), o₁:o₂)...)
    append!(v, 10. ^ (o₂ + 1))
    return v
end

##-----------------------------------------------------------------------------
# STAGE 1

ens = ensemble(
    product(τ, σ),
    t₁,
    t₂,
    nwide,
    nstep,
    nstore,
    𝒻W
)

saveensemble(
    datadir(
        "sims",
        "ensemble_wide.jld2"
    ),
    ens
)

#find where temperature range is appreciable and peak fCO2 isn't wacky high
df = filter(
    [:Trange,:fmax] => (a,b) -> (a > 5) & (b < 9e5),
    combine(
        groupby(
            framevariable(:T, ens),
            [:τ, :σ]
        ),
        [:Tmax,:Tmin] => ((a,b) -> median(a - b)) => :Trange,
        :fmax => median => :fmax
    )
)

##-----------------------------------------------------------------------------
# STAGE 2

#simulate ensemble and save directly
saveensemble(
    datadir(
        "sims",
        "ensemble_deep.jld2"
    ),
    ensemble(
        zip(df.τ, df.σ),
        t₁,
        t₂,
        ndeep,
        nstep,
        nstore,
        𝒻W
    )
)
