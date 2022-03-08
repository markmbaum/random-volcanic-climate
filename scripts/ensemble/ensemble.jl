using DrWatson
@quickactivate "Random Volcanic Climate"
using Pkg
Pkg.instantiate()

push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using IterTools: product
using AxisArrays
using MultiAssign
using Base.Threads: @threads, nthreads, threadid
using ProgressMeter

## FUNCTIONS

function ensemble(params,
                  t₁::Float64,
                  t₂::Float64,
                  nrealize::Int,
                  nstep::Int,
                  nstore::Int,
                  𝒻W::F
                  ) where {F<:Function}
    println(stdout, "starting ensemble with $(nthreads()) threads")
    #number of parameter combinations
    L = length(params)
    #total number of simulations
    N = L*nrealize
    println(stdout, "$N total simulations")
    flush(stdout)
    #predict the time samples and their indices
    idx = Int.(round.(range(1, nstep, nstore)))
    tₛ = round.(LinRange(t₁, t₂, nstep)[idx], sigdigits=4)
    #allocate an array for the parameter combinations and fill in values
    p = AxisArray(zeros(Float32, 2, N), parameter=[:τ, :σ], trial=1:N)
    i = 1
    for (τ, σ) ∈ params, _ ∈ 1:nrealize
        p[:,i] .= τ, σ
        i += 1
    end
    #allocate an array for carbon and outgassing at all stored times
    res = AxisArray(
        zeros(Float32, nstore, N, 2),
        time=tₛ,
        trial=1:N,
        var=[:C, :V]
    )
    #space for all steps of in-place simulations
    @multiassign c, v = zeros(nstep, nthreads())
    #initial carbon reservoir size
    C₁ = 𝒻Cₑ(t₁)
    #initial outgassing rate, subject to spinup
    V₁ = Vᵣ
    #simulate
    progress = Progress(N)
    @threads for i ∈ 1:N
        id = threadid()
        simulate!(
            view(c, :, id),
            view(v, :, id),
            t₁,
            t₂,
            C₁,
            V₁,
            𝒻W,
            initparams(
                τ=p[1,i],
                σ=p[2,i]
            )
        )
        #store selected values
        res[:,i,:C] .= @view c[idx,id]
        res[:,i,:V] .= @view v[idx,id]
        #progress updates
        next!(progress)
    end
    return p, res
end

## INPUT/PARAMETERS

#simulation start time [Gyr]
t₁ = 2.5
#simulation end time [Gyr]
t₂ = 4.5
#values for outgassing relaxation
τ = exp10.(LinRange(6, 9, 4))
#values for outgassing variance
σ = exp10.(LinRange(-6, -4, 4))
#weathering function
𝒻W(C,t) = 𝒻whak(C, t, β=0)
#number of simulations per parameter combination
nrealize = 20*nthreads()
#number of steps for each simulation
nstep = 1_000_000
#number of time slices to store
nstore = 51

## MAIN

#create parameter combinations
params = product(τ, σ)

#simulate
p, res = ensemble(
    params,
    Float64(t₁),
    Float64(t₂),
    nrealize,
    nstep,
    nstore,
    𝒻W
)

##

safesave(
    datadir(
        "sims",
        "ensemble.jld2"
    ),
    Dict(
        "p"=>p,
        "res"=>res
    )
)
