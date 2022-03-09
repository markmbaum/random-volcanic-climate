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
    t = round.(LinRange(t₁, t₂, nstep+1)[idx], sigdigits=4)
    #allocate arrays for the parameter combinations and fill in values
    @multiassign τ, σ = zeros(N)
    i = 1
    for p ∈ params, _ ∈ 1:nrealize
        τ[i] = p[1]
        σ[i] = p[2]
        i += 1
    end
    #allocate an array for carbon and outgassing at all stored times
    res = AxisArray(
        zeros(Float32, 4, nstore, N),
        var=[:C, :V, :T, :W],
        time=t,
        trial=1:N
    )
    #space for all steps of in-place simulations
    @multiassign c, v = zeros(nstep, nthreads())
    #initial carbon reservoir size
    C₁ = 𝒻Cₑ(t₁)
    #initial outgassing rate, subject to spinup
    V₁ = Vᵣ
    #simulate
    progress = Progress(N, output=stdout)
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
                τ=τ[i],
                σ=σ[i]
            )
        )
        #store selected values
        res[:C,:,i] .= @view c[idx,id]
        res[:V,:,i] .= @view v[idx,id]
        #also store temperature and weathering
        res[:T,:,i] .= C2T.(view(res,:C,:,i), t)
        res[:W,:,i] .= 𝒻W.(view(res,:C,:,i), t)
        #progress updates
        next!(progress)
    end
    return t, τ, σ, res
end

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

## MAIN

#create parameter combinations
params = product(τ, σ)

#simulate
t, τ, σ, res = ensemble(
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
        "t"=>t,
        "τ"=>τ,
        "σ"=>σ,
        "res"=>res
    )
)
