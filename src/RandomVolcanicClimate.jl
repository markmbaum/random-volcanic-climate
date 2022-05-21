module RandomVolcanicClimate

using Base.Threads: @threads, nthreads, threadid
using Roots
using BasicInterpolators: ChebyshevInterpolator
using ForwardDiff: derivative
using GEOCLIM: godderis, whak, mac
using UnPack
using MultiAssign
using AxisArrays
using DrWatson
using DataFrames: DataFrame, insertcols!, stack
using Statistics: median
using ProgressMeter: Progress, next!

#------------------------------------------------------------------------------
# time units/conversions

export day, yr, Myr, Gyr

const day = 60*60*24 #seconds in a day
const yr = 365.25*day
const Myr = 1e6*yr
const Gyr = 1e3*Myr

#------------------------------------------------------------------------------
# immutable physical constants

export 𝐭, 𝐠, 𝛍, 𝐑ₑ, 𝐒ₑ, 𝐅

#age of solar system [Gyr]
const 𝐭 = 4.5
#surface gravity [m/s^2]
const 𝐠 = 9.8
#molar mass of CO2 [kg/mole]
const 𝛍 = 0.044
#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6
#the Earth's surface area [m^2]
const 𝐒ₑ = 4π*𝐑ₑ^2
#solar "constant" [W/m^2]
const 𝐅 = 1366.0

#------------------------------------------------------------------------------
# physical constants of choice, reference values

export fCO2ᵣ, Cᵣ, Tᵣ, pᵣ, αᵣ, OLRᵣ, hᵣ, pᵣ, Vᵣ, γ

#reference molar concentration of CO2 [ppm]
const fCO2ᵣ = 285.0
#reference total carbon in ocean-atmosphere [teramole]
const Cᵣ = 3.193e18/1e12
#reference temperature [K]
const Tᵣ = 288.0
#reference precipitation rate [m/s]
const pᵣ = 1.0/yr
#proportionality between precipitation and runoff
const Γ = 0.2
#fractional change in precipitation per K change in temperature
const ϵ = 0.03
#reference albedo [-]
const αᵣ = 0.3
#reference OLR [W/m^2]
const OLRᵣ = (1 - αᵣ)*𝐅/4
#default parameter for CO2 ocean-atmosphere partioning [teramole]
const hᵣ = 2.3269250670587494e20/1e12
#reference atmospheric pressure excluding CO2
const Pᵣ = 1e5 - 28.5
#default volcanic CO2 outgassing rate [teramole/yr]
const Vᵣ = 7.0
#land fraction for weathering [-]
const γ = 0.3

#OLR response to temperature
const a = 2.0
#OLR response to pCO2
const b = 5.35

#------------------------------------------------------------------------------
# component physical equations

export 𝒻☉, 𝒻F, 𝒻S, 𝒻T, 𝒻ϕ, 𝒻pCO2, 𝒻fCO2, C2T, 𝒻OLR, 𝒻RI, 𝒻L, 𝒻p, 𝒻q

#stellar luminosity fraction over time [Gya]
𝒻☉(t=𝐭) = 1/(1 + (2/5)*(1 - t/𝐭))

#instellation over time [W/m^2]
𝒻F(t=𝐭) = 𝒻☉(t)*𝐅

#surface area averaged radiation [W/m^2]
𝒻S(t=𝐭) = 𝒻F(t)/4

#global temperature [K]
𝒻T(fCO2=fCO2ᵣ, t=𝐭, α=αᵣ) = ((1 - α)*𝒻S(t) - OLRᵣ + b*log(fCO2/fCO2ᵣ))/a + Tᵣ

#fraction of carbon in the atmosphere [-]
# Mills, Benjamin, et al. "Timing of Neoproterozoic glaciations linked to transport-limited global weathering." Nature geoscience 4.12 (2011): 861-864.
𝒻ϕ(C=Cᵣ, h=hᵣ) = 0.78*C/(C + h)

#partial pressure of CO2 [Pa]
𝒻pCO2(C=Cᵣ, h=hᵣ) = 𝒻ϕ(C,h)*(C*1e12)*𝛍*𝐠/𝐒ₑ

#molar concentration of CO2 [ppmv]
function 𝒻fCO2(C=Cᵣ, P=Pᵣ, h=hᵣ)
    @assert C > 0
    pCO2 = 𝒻pCO2(C, h)
    return 1e6*pCO2/(pCO2 + P)
end

#convert carbon reservior [Tmole] directly to temperature [K]
C2T(C, t=𝐭, α=αᵣ) = 𝒻T(𝒻fCO2(C), t, α)

#outgoing longwave radiation [W/m^2]
𝒻OLR(T=Tᵣ, fCO2=fCO2ᵣ) = OLRᵣ + a*(T - Tᵣ) - b*log(fCO2/fCO2ᵣ)

#radiative imbalance [W/m^2]
𝒻RI(T=Tᵣ, fCO2=fCO2ᵣ, t=𝐭, α=αᵣ) = (1 - α)*𝒻S(t) - 𝒻OLR(T, fCO2)

#latent heat of vaporization of water [J/m^3]
𝒻L(T=Tᵣ) = 1.918e9*(T/(T - 33.91))^2

#global precipitation [m/s]
𝒻p(T=Tᵣ, t=𝐭, α=αᵣ) = max(min(pᵣ*(1 + ϵ*(T - Tᵣ)), (1 - α)*𝒻S(t)/𝒻L(T)), zero(T))

#global runoff [m/s]
𝒻q(T=Tᵣ, t=𝐭, α=αᵣ) = Γ*𝒻p(T, t, α)

#------------------------------------------------------------------------------
# equilibrium carbon content over time and its derivative

export 𝒻Cₑ, Χ, dΧ

function 𝒻Cₑ(t=𝐭, T=Tᵣ)
    #find the root in log space because total carbon is a big number
    Roots.find_zero(
        x -> 𝒻RI(T, 𝒻fCO2(exp10(x)), t),
        (-10, 10) #bracketing initial guesses in log space
    ) |> exp10
end

#simple struct to rapidly interpolate χ values instead of root finding each time
struct Χ #capital Chi here
    interpolator::ChebyshevInterpolator{32,Float64}
end

#constructor
function Χ(Tₑ::Real=Tᵣ) 
    ChebyshevInterpolator(
        t -> log(𝒻Cₑ(t, Float64(Tₑ))), #function to approximate
        2.5, #time interval beginning
        𝐭, #time interval end
        32 #number of interpolation nodes, 32 is more than enough
    ) |> Χ
end

#functor for carbon reservior yielding equilibrium temperature
(χ::Χ)(t) = exp(χ.interpolator(t))

#same idea, but for rate of carbon change to maintain equilibrium
struct dΧ #capital Chi here
    interpolator::ChebyshevInterpolator{32,Float64}
end

function dΧ(Tₑ::Real=Tᵣ)
    #construct the carbon curve before taking its derivative
    χ = Χ(Tₑ)
    #make an interpolator for its derivative using ForwardDiff.jl
    ChebyshevInterpolator(
        t -> log(-derivative(χ, t)), #function to approximate
        2.5, #time interval beginning
        𝐭, #time interval end
        32 #number of interpolation nodes, 32 is more than enough
    ) |> dΧ
end

(dχ::dΧ)(t) = -exp(dχ.interpolator(t))

#------------------------------------------------------------------------------
# weathering and some equilibrium finders

export 𝒻whak, 𝒻mac, 𝒻Wₑ, 𝒻Tₑ

function preweathering(C, t)
    fCO2 = 𝒻fCO2(C) #CO2 concentration [ppm]
    T = 𝒻T(fCO2, t) #global temperature [K]
    q = Γ*pᵣ #global runoff [m/s]
    return fCO2, T, q
end

function 𝒻whak(C=Cᵣ, t=𝐭; k=0.2287292550091995, β=0.0)
    fCO2, T, q = preweathering(C, t)
    #weathering rate [mole/second/m²]
    w = whak(q, T, fCO2, k, 11.1, Tᵣ, fCO2ᵣ, β)
    #global weathering [teramole/year]
    w*(0.3*𝐒ₑ*yr/1e12)
end

function 𝒻mac(C=Cᵣ, t=𝐭; Λ=6.1837709746872e-5, β=0.2)
    fCO2, T, q = preweathering(C, t)
    #weathering rate [mole/second/m²]
    w = mac(q, T, fCO2, 11.1, Tᵣ, fCO2ᵣ, Λ=Λ, β=β)
    #global weathering [teramole/year]
    w*(0.3*𝐒ₑ*yr/1e12)
end

#finds carbon reservoir [teramole] where weathering balances volcanism [teramole/yr]
function 𝒻Wₑ(𝒻W::F, V=Vᵣ, t=𝐭) where {F<:Function}
    Roots.find_zero(
        x -> 𝒻W(exp10(x),t) - V,
        (-10, 15),
        Roots.Brent()
    ) |> exp10
end

#finds temperature [K] where weathering balances volcanism [teramole/yr]
#this can be done analytically for simple weathering formulae (whak with β=0)
𝒻Tₑ(𝒻W::F, V=Vᵣ, t=𝐭) where {F<:Function} = 𝒻T(𝒻Wₑ(𝒻W, V, t) |> 𝒻fCO2, t)

#------------------------------------------------------------------------------
# integration/modeling

export initparams
export simulate!, simulate

function initparams(;
    μ::Real=Vᵣ, #mean volcanic outgassing rate [teramole/yr]
    τ::Real=1e7, #outgassing relaxation timescale [yr]
    σ::Real=1e-4, #outgassing std [teramole/yr]
    Vₘ::Real=0.0, #minimum outgassing rate [teramole/yr]
    Cₘ::Real=Cᵣ/1e6, #minimum allowable C reservoir size [teramole]
    spinup::Real=0.5 #spinup time [Gyr]
    )::NamedTuple
    (
        μ=Float64(μ),
        τ=Float64(τ),
        σ=Float64(σ),
        Vₘ=Float64(Vₘ),
        Cₘ=Float64(Cₘ),
        spinup=Float64(spinup)
    )
end

function step(tᵢ, Cᵢ, Vᵢ, Δt, 𝒻W::F, params) where {F<:Function}
    @unpack μ, τ, σ, Vₘ, Cₘ = params
    Cᵢ₊₁ = Cᵢ + Δt*(Vᵢ - 𝒻W(Cᵢ, tᵢ))
    Vᵢ₊₁ = Vᵢ + Δt*(μ - Vᵢ)/τ + √(Δt)*σ*randn()
    return max(Cᵢ₊₁, Cₘ), max(Vᵢ₊₁, Vₘ)
end

function simulate!(C::AbstractVector,
                   V::AbstractVector,
                   t₁::Float64,
                   t₂::Float64,
                   C₁::Float64,
                   V₁::Float64,
                   𝒻W::F,
                   params::NamedTuple
                   )::Nothing where {F<:Function}
    #check output lengths
    @assert length(C) == length(V)
    nstep = length(C) - 1
    #initialize time stepping
    t = LinRange(t₁, t₂, nstep+1) #units of Gyr
    Δt = 1e9*(t₂ - t₁)/nstep #units of yr
    #spin up
    spinup = 1e9*params[:spinup] #convert to Gyr
    tspin = 0.0
    while tspin < spinup
        C₁, V₁ = step(t₁, C₁, V₁, Δt, 𝒻W, params)
        tspin += Δt
    end
    #initial values
    C[1] = C₁
    V[1] = V₁
    #solve/integrate
    @inbounds for i ∈ 2:nstep+1
        C[i], V[i] = step(t[i-1], C[i-1], V[i-1], Δt, 𝒻W, params)
    end
    return nothing
end

function simulate(params=initparams()::NamedTuple;
                  t₁=2.5,
                  t₂=4.5,
                  C₁=nothing,
                  V₁=nothing,
                  𝒻W=𝒻whak,
                  nstep=1_000_000)
    C = zeros(nstep + 1)
    V = zeros(nstep + 1)
    t = LinRange(t₁, t₂, nstep + 1)
    simulate!(
        C,
        V,
        Float64(t₁),
        Float64(t₂),
        Float64(isnothing(C₁) ? 𝒻Cₑ(t₁) : C₁),
        Float64(isnothing(V₁) ? Vᵣ : V₁),
        𝒻W,
        params
    )
    return t, C, V
end

#------------------------------------------------------------------------------
# main ensemble function

export ensemble, snowballtime

function snowballtime(t::AbstractVector,
                      T::AbstractVector,
                      Tsnow::Real=280f0
                      )::Float32
    @assert length(t) == length(T)
    if any(x -> x < Tsnow, T)
        i = findfirst(x -> x < Tsnow, T)
        if i == 1
            return t[1]
        else
            #linear interpolation
            @inbounds begin
                t₁ = t[i-1]
                t₂ = t[i]
                T₁ = T[i-1]
                T₂ = T[i]
            end
            return (Tsnow - T₁)*(t₂ - t₁)/(T₂ - T₁) + t₁
        end
    end
    return NaN32
end

#barrier
function ensemble(params, t₁, t₂, nrealize, nstep, nstore, 𝒻W::F) where {F<:Function}
    ensemble(
        params,
        Float64(t₁),
        Float64(t₂),
        Int64(nrealize),
        Int64(nstep),
        Int64(nstore),
        𝒻W
    )
end

#params should be an iterable of (τ, σ) containers
function ensemble(params,
                  t₁::Float64,
                  t₂::Float64,
                  nrealize::Int,
                  nstep::Int,
                  nstore::Int,
                  𝒻W::F,
                  Tsnow::Float64=280.0,
                  ) where {F<:Function}
    println(stdout, "starting ensemble with $(nthreads()) threads")
    #number of parameter combinations
    L = length(params)
    #total number of simulations
    N = L*nrealize
    println(stdout, "$N total simulations")
    flush(stdout)
    #time samples along the dense time steps of each simulation
    tsim = LinRange(t₁, t₂, nstep + 1)
    #time samples and their indices
    idx = Int.(round.(range(1, nstep, nstore)))
    tstore = round.(tsim[idx], sigdigits=4)
    #allocate arrays for the parameter combinations and fill in values
    @multiassign τ, σ = zeros(Float32, N)
    i = 1
    for p ∈ params, _ ∈ 1:nrealize
        τ[i] = p[1]
        σ[i] = p[2]
        i += 1
    end
    #allocate an array for carbon, outgassing, and prognostics at stored times
    R = AxisArray(
        zeros(Float32, 4, nstore, N),
        var=[:C, :V, :T, :W],
        time=tstore,
        trial=1:N
    )
    #space for columns representing a single result for each trial
    S = AxisArray(
        zeros(Float32, 7, N),
        var=[:tsno, :Cmax, :Vmax, :Tmax, :Tmin, :tmax, :tmin],
        trial=1:N
    )
    #space for all steps of in-place simulations
    @multiassign c, v = zeros(nstep + 1, nthreads())
    #initial carbon reservoir size
    C₁ = 𝒻Cₑ(t₁)
    #initial outgassing rate, subject to spinup
    V₁ = Vᵣ
    #simulate
    progress = Progress(N, output=stdout)
    @threads for i ∈ 1:N
        id = threadid()
        cᵢ = @view c[:,id]
        vᵢ = @view v[:,id]
        simulate!(cᵢ, vᵢ, t₁, t₂, C₁, V₁, 𝒻W,
            initparams(
                τ=τ[i],
                σ=σ[i]
            )
        )
        #store main variables at time slices
        R[:C,:,i] .= @view c[idx,id]
        R[:V,:,i] .= @view v[idx,id]
        #also store temperature and weathering
        R[:T,:,i] .= C2T.(view(R,:C,:,i), tstore)
        R[:W,:,i] .= 𝒻W.(view(R,:C,:,i), tstore)
        #full temperature solution
        Tsim = C2T.(cᵢ, tsim)
        #time to snowball
        S[:tsno,i] = snowballtime(tsim, Tsim, Tsnow)
        #extreme C, V, and T values
        S[:Cmax,i] = maximum(cᵢ)
        S[:Vmax,i] = maximum(vᵢ)
        S[:Tmax,i], j = findmax(Tsim)
        S[:tmax,i] = tsim[j]
        S[:Tmin,i], j = findmin(Tsim)
        S[:tmin,i] = tsim[j]
        #progress updates
        next!(progress)
    end
    return tstore, τ, σ, R, S
end

#------------------------------------------------------------------------------
# some handy functions for saving, loading, organizing ensemble results

export saveensemble, loadensemble, framevariable, frameensemble, timecols, stacktimes

function saveensemble(fn, t, τ, σ, R, S)
    safesave(fn, Dict("t"=>t, "τ"=>τ, "σ"=>σ, "R"=>R, "S"=>S))
end

saveensemble(fn, X) = saveensemble(fn, X...)

function loadensemble(fn::String)
    ens = wload(fn)
    @unpack t, τ, σ, R, S = ens
    return (t=t, τ=τ, σ=σ, R=R, S=S)
end

function framevariable(X::AxisArray, t, τ, σ, S)
    N = size(X, 2)
    L = length(t)
    cols = [:τ, :σ, :fmax, :tsno, :Cmax, :Vmax, :Tmax, :Tmin, :tmax, :tmin]
    iₜ = length(cols)
    df = DataFrame(
        zeros(Float32, N, length(t) + iₜ),
        vcat(
            cols,
            1:L .|> Symbol
        )
    )
    #main parameters
    df[:,:τ] = τ
    df[:,:σ] = σ
    #derived fCO2 peak
    df[:,:fmax] = 𝒻fCO2.(S[:Cmax,:])
    #all the singles
    for col ∈ cols[4:end]
        df[:,col] = S[col,:]
    end
    #snapshots of main variable
    df[:,iₜ+1:end] = X[:,:]'
    return df
end

function framevariable(X::AxisArray, ens::NamedTuple)
    @unpack t, τ, σ, S = ens
    framevariable(X, t, τ, σ, S)
end

framevariable(var::Symbol, t, τ, σ, R, S) = framevariable(R[var], t, τ, σ, S)

framevariable(var::Symbol, X) = framevariable(var, X...)

function frameensemble(t, τ, σ, R, S)
    var = (:C, :V, :T, :W)
    dfs = (framevariable(v, t, τ, σ, R, S) for v ∈ var)
    return t, (; zip(var, dfs)...)
end

frameensemble(X) = frameensemble(X...)

timecols(df::DataFrame) = [x for x ∈ names(df) if !isnothing(tryparse(Int, x))]

function stacktimes(df)
    #stack/melt all the time columns
    stack(
        df,
        df |> timecols .|> Symbol,
        variable_name="time"
    )
end

function 𝒻Tₑ(𝒻W::F, V::AxisArray) where {F<:Function}
    #time slices
    t = V.axes[1] |> collect
    #new array
    Tₑ = similar(V)
    N, M = size(V)
    #fill it in
    p = Progress(M)
    @threads for j ∈ 1:M
        @inbounds for i ∈ 1:N
            try
                Tₑ[i,j] = 𝒻Tₑ(𝒻W, V[i,j], t[i])
            catch err
                Tₑ[i,j] = NaN
            end
        end
        next!(p)
    end
    return Tₑ
end

#------------------------------------------------------------------------------
# misc

export 𝒻gya

𝒻gya(t) = 𝐭 - t

end