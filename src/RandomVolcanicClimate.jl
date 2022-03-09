module RandomVolcanicClimate

using Base.Threads: @threads
using Roots: find_zero, Newton
using BasicInterpolators: ChebyshevInterpolator
using ForwardDiff: derivative
using GEOCLIM: godderis, whak, mac
using UnPack
using MultiAssign
using AxisArrays
using DrWatson

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
#default parameter for CO2 partioning [teramole]
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

#absorbed radiation [W/m^2]
𝒻S(t=𝐭, α=αᵣ) = (1 - α)*𝒻F(t)/4

#global temperature [K]
𝒻T(fCO2=fCO2ᵣ, t=𝐭, α=αᵣ) = (𝒻S(t,α) - OLRᵣ + b*log(fCO2/fCO2ᵣ))/a + Tᵣ

#fraction of carbon in the atmosphere [-]
# Mills, Benjamin, et al. "Timing of Neoproterozoic glaciations linked to transport-limited global weathering." Nature geoscience 4.12 (2011): 861-864.
𝒻ϕ(C=Cᵣ, h=hᵣ) = 0.78*C/(C + h)

#partial pressure of CO2 [Pa]
𝒻pCO2(C=Cᵣ, h=hᵣ) = 𝒻ϕ(C,h)*(C*1e12)*𝛍*𝐠/𝐒ₑ

#molar concentration of CO2 [ppmv]
function 𝒻fCO2(C=Cᵣ, P=Pᵣ, h=hᵣ)
    pCO2 = 𝒻pCO2(C, h)
    return 1e6*pCO2/(pCO2 + P)
end

#convert carbon reservior [Tmole] directly to temperature [K]
function C2T(C, t=𝐭, α=αᵣ)
    @assert C > 0
    𝒻T(𝒻fCO2(C), t, α)
end

#outgoing longwave radiation [W/m^2]
𝒻OLR(T=Tᵣ, fCO2=fCO2ᵣ) = OLRᵣ + a*(T - Tᵣ) - b*log(fCO2/fCO2ᵣ)

#radiative imbalance [W/m^2]
𝒻RI(T=Tᵣ, fCO2=fCO2ᵣ, t=𝐭, α=αᵣ) = 𝒻S(t, α) - 𝒻OLR(T, fCO2)

#latent heat of vaporization of water [J/m^3]
𝒻L(T=Tᵣ) = 1.918e9*(T/(T - 33.91))^2

#global precipitation [m/s]
𝒻p(T=Tᵣ, t=𝐭) = min(pᵣ*(1 + ϵ*(T - Tᵣ)), 𝒻S(t)/𝒻L(T))

#global runoff [m/s]
𝒻q(T=Tᵣ, t=𝐭) = Γ*𝒻p(T,t)

#------------------------------------------------------------------------------
# equilibrium carbon content over time and its derivative

export 𝒻Cₑ, Χ, dΧ

function 𝒻Cₑ(t=𝐭, T=Tᵣ)
    #find the root in log space because total carbon is a big number
    T = exp10(
        find_zero(
            x -> 𝒻RI(T, 𝒻fCO2(exp10(x)), t),
            (-10, 10) #good bracketing initial guesses in log space
        )
    )
    return T
end

#simple struct to rapidly interpolate χ values instead of root finding each time
struct Χ #capital Chi here
    interpolator::ChebyshevInterpolator{32,Float64}
end

#constructor
function Χ(Tₑ::Real=Tᵣ) 
    I = ChebyshevInterpolator(
        t -> log(𝒻Cₑ(t, Float64(Tₑ))), #function to approximate
        2.5, #time interval beginning
        𝐭, #time interval end
        32 #number of interpolation nodes, 32 is more than enough
    )
    Χ(I)
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
    I = ChebyshevInterpolator(
        t -> log(-derivative(χ, t)), #function to approximate
        2.5, #time interval beginning
        𝐭, #time interval end
        32 #number of interpolation nodes, 32 is more than enough
    )
    dΧ(I)
end

(dχ::dΧ)(t) = -exp(dχ.interpolator(t))

#------------------------------------------------------------------------------
# weathering

export 𝒻whak, 𝒻mac, 𝒻Wₑ

function preweathering(C, t)
    fCO2 = 𝒻fCO2(C) #CO2 concentration [ppm]
    T = 𝒻T(fCO2, t) #global temperature [K]
    q = 𝒻q(Tᵣ, t) #global runoff [m/s]
    return fCO2, T, q
end

function 𝒻whak(C=Cᵣ, t=𝐭; k=0.2287292550091995, β=0.0)
    fCO2, T, q = preweathering(C, t)
    #weathering rate [mole/second]
    w = whak(q, T, fCO2, k, 11.1, Tᵣ, fCO2ᵣ, β)
    #global weathering [teramole/year]
    w*(0.3*𝐒ₑ*yr/1e12)
end

function 𝒻mac(C=Cᵣ, t=𝐭; Λ=6.1837709746872e-5, β=0.2)
    fCO2, T, q = preweathering(C, t)
    #weathering rate [mole/second]
    w = mac(q, T, fCO2, 11.1, Tᵣ, fCO2ᵣ, Λ=Λ, β=β)
    #global weathering [teramole/year]
    w*(0.3*𝐒ₑ*yr/1e12)
end

#finds carbon reservoir where weathering balances volcanism
function 𝒻Wₑ(𝒻W::F, V=Vᵣ, t=𝐭) where {F}
    C = exp10(
        find_zero(
            x -> 𝒻W(exp10(x),t) - V,
            (-10, 10)
        )
    )
    return C
end

#------------------------------------------------------------------------------
# integration/modeling

export initparams
export simulate!, simulate

function initparams(;
    μ::Real=Vᵣ, #mean volcanic outgassing rate [teramole/yr]
    τ::Real=1e7, #outgassing relaxation timescale [yr]
    σ::Real=1e-4, #outgassing variance []
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
    Cᵢ₊₁ = max( Cᵢ + Δt*(Vᵢ - 𝒻W(Cᵢ, tᵢ)), Cₘ )
    Vᵢ₊₁ = max( Vᵢ + Δt*(μ - Vᵢ)/τ + √(Δt)*σ*randn(), Vₘ)
    return Cᵢ₊₁, Vᵢ₊₁
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
    t = LinRange(t₁, t₂, nstep+1)
    Δt = 1e9*(t₂ - t₁)/nstep
    #spin up
    @unpack spinup = params
    spinup *= 1e9
    tspin = 0.0
    while tspin < spinup
        C₁, V₁ = step(t₁, C₁, V₁, Δt, 𝒻W, params)
        tspin += Δt
    end
    #initial values
    C[1] = C₁
    V[1] = V₁
    #solve/integrate
    for i ∈ 2:nstep+1
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

end