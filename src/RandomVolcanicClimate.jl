module RandomVolcanicClimate

using Base.Threads: @threads
using Distributions
using Roots: find_zero
using BasicInterpolators: ChebyshevInterpolator
using ForwardDiff: derivative
using GEOCLIM: whak, mac

include("PowerLawDistribution.jl")
using .PowerLawDistribution
export PowerLaw, ∫, mean, var

#------------------------------------------------------------------------------
# time units/conversions

export day, yr, Myr, Gyr

const day = 60*60*24 #seconds in a day
const yr = 365.25*day
const Myr = 1e6*yr
const Gyr = 1e3*Myr

#------------------------------------------------------------------------------
# immutable physical constants

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

export fCO2ᵣ, Cᵣ, Tᵣ, pᵣ, αᵣ, OLRᵣ, hᵣ, pᵣ

#reference molar concentration of CO2 [ppm]
const fCO2ᵣ = 285.0
#reference total carbon in ocean-atmosphere [moles]
const Cᵣ = 3.193e18
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
#default parameter for CO2 partioning
const hᵣ = 2.3269250670587494e20
#reference atmospheric pressure excluding CO2
const Pᵣ = 1e5 - 28.5

#OLR response to temperature
const a = 2.0
#OLR response to pCO2
const b = 5.35

#------------------------------------------------------------------------------
# component physical equations

export 𝒻☉, 𝒻F, 𝒻S, 𝒻T, 𝒻ϕ, 𝒻pCO2, 𝒻fCO2, 𝒻OLR, 𝒻RI, 𝒻L, 𝒻p, 𝒻q

#stellar luminosity fraction over time [Gya]
𝒻☉(t=𝐭) = 1/(1 + (2/5)*(1 - t/𝐭))

#instellation over time [W/m^2]
𝒻F(t=𝐭) = 𝒻☉(t)*𝐅

#absorbed radiation [W/m^2]
𝒻S(t=𝐭, α=αᵣ) = (1 - α)*𝒻F(t)/4

#global temperature [K]
𝒻T(fCO2=fCO2ᵣ, t=𝐭, α=αᵣ) = (𝒻S(t,α) - OLRᵣ + b*log(fCO2/fCO2ᵣ))/a + Tᵣ

#fraction of carbon in the atmosphere [-], see Mills et al., 2011
𝒻ϕ(C=Cᵣ, h=hᵣ) = 0.78*C/(C + h)

#partial pressure of CO2 [Pa]
𝒻pCO2(C=Cᵣ, h=hᵣ) = 𝒻ϕ(C,h)*C*𝛍*𝐠/𝐒ₑ

#molar concentration of CO2 [ppmv]
function 𝒻fCO2(C=Cᵣ, P=Pᵣ, h=hᵣ)
    pCO2 = 𝒻pCO2(C, h)
    return 1e6*pCO2/(pCO2 + P)
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

function 𝒻Cₑ(t=𝐭, Tₑ=Tᵣ)
    #find the root in log space because total carbon is a big number
    exp10(
        find_zero(
            x -> 𝒻RI(Tₑ, 𝒻fCO2(exp10(x)), t),
            (16, 22) #good bracketing initial guesses
        )
    )
end

#simple struct to rapidly interpolate χ values instead of root finding
struct Χ{𝒯<:ChebyshevInterpolator} #capital Chi here
    interpolator::𝒯
end

function Χ(Tₑ=Tᵣ; N::Int=32)
    I = ChebyshevInterpolator(
        t -> log(𝒻Cₑ(t, Float64(Tₑ))), #function to approximate
        2.5, #time interval beginning
        𝐭, #time interval end
        N #number of interpolation nodes
    )
    Χ(I)
end

#functor for carbon reservior yielding equilibrium temperature
(χ::Χ)(t) = exp(χ.interpolator(t))

#same idea, but for rate of carbon change to maintain equilibrium
struct dΧ{𝒯<:ChebyshevInterpolator} #capital Chi here
    interpolator::𝒯
end

function dΧ(Tₑ=Tᵣ; N::Int=32)
    #construct the carbon curve
    χ = Χ(Tₑ, N=N)
    #make an interpolator for its derivative using ForwardDiff.jl
    I = ChebyshevInterpolator(
        t -> log(-derivative(χ, t)), #function to approximate
        2.5, #time interval beginning
        𝐭, #time interval end
        N #number of interpolation nodes
    )
    dΧ(I)
end

(dχ::dΧ)(t) = -exp(dχ.interpolator(t))

#------------------------------------------------------------------------------
# weathering


#------------------------------------------------------------------------------
# integration/modeling

export step
export integrate, integrations
export simulate, simulations

function setup(V, t₁, t₂, nstep)
    @assert (t₁ > 0) & (t₂ > 0)
    @assert t₂ > t₁
    @assert nstep > 0
    t = t₁
    Δt = (t₂ - t₁)/nstep
    Δtₛ = √(Δt)
    μ = mean(V)
    return t, Δt, Δtₛ, μ
end

function step(t, C, Δt, Δtₛ, μ, V, 𝒻W)::Float64
    #ordinary part
    C += Δt*(μ - 𝒻W(t,C))
    #random part
    C += Δtₛ*(rand(V) - μ)*8e5
    return C
end

#type restricted function
function simulate(V::Sampleable{Univariate,Continuous},
                  𝒻W::F,
                  t₁::Float64,
                  t₂::Float64,
                  C₁::Float64,
                  nstep::Int) where {F<:Function}
    t, Δt, Δtₛ, μ = setup(V, t₁, t₂, nstep)
    t = LinRange(t₁, t₂, nstep+1)
    C = zeros(nstep+1)
    C[1] = C₁
    for i ∈ 2:nstep+1
        C[i] = step(t[i-1], C[i-1], Δt, Δtₛ, μ, V, 𝒻W)
    end
    return t, C
end

#function barrier, type flexible
function simulate(V, 𝒻W, t₁=2.5, t₂=4.5; nstep::Int=100_000)
    simulate(V, 𝒻W, Float64(t₁), Float64(t₂), Float64(𝒻Cₑ(t₁)), nstep)
end

function integrate(V::Sampleable{Univariate,Continuous},
                   𝒻W::F,
                   C₁::Float64,
                   t₁::Float64,
                   t₂::Float64,
                   nstep::Int) where {F<:Function}
    t, Δt, Δtₛ, μ = setup(V, t₁, t₂, nstep)
    C = C₁
    for _ ∈ 1:nstep
        C = step(t, C, Δt, Δtₛ, μ, V, 𝒻W)
        t += Δt
    end
    return C
end

function integrate(V, 𝒻W, t₁=2.5, t₂=4.5; nstep::Int=10_000)
    integrate(V, 𝒻W, Float64(t₁), Float64(t₂), 𝒻Cₑ(t₁), nstep)
end

function integrations(N::Int, args...; kw...)
    C = zeros(N)
    @threads for i ∈ 1:N
        C[i] = integrate(args...; kw...)
    end
    return C
end

end