module RandomVolcanicClimate

using Roots
using UnPack
using Random: AbstractRNG
using Cubature: hquadrature
using Base.Threads: @threads

include("power_law.jl")

#------------------------------------------------------------------------------
# physical constants

#age of solar system
const 𝐭 = 4.5
#surface gravity [m/s^2]
const 𝐠 = 9.8
#reference instellation (solar "constant") [W/m^2]
const F₀ = 1366.0
#reference molar concentration of CO2 [ppm]
const fCO2₀ = 280
#reference temperature [K]
const T₀ = 288.0
#reference albedo [-]
const α₀ = 0.3
#reference OLR [W/m^2]
const OLR₀ = (1 - α₀)*F₀/4
#OLR response to pCO2
const b = 5.35
#default parameter for CO2 partioning
const 𝐡 = 2.3269250670587494e20
#molar mass of CO2 [kg/mole]
const 𝛍 = 0.044
#reference atmospheric pressure excluding CO2
const 𝐏₀ = 1e5 - 28.5
#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6
#the Earth's surface area [m^2]
const 𝐒ₑ = 4π*𝐑ₑ^2

#------------------------------------------------------------------------------
# component physical equations

export 𝒻☉, 𝒻F, 𝒻S, 𝒻T, 𝒻ϕ, 𝒻pCO2, 𝒻fCO2

#stellar luminosity fraction over time [Gya]
𝒻☉(t) = 1/(1 + (2/5)*(1 - t/𝐭))

#instellation over time [W/m^2]
𝒻F(t) = 𝒻☉(t)*F₀

#absorbed radiation [W/m^2]
𝒻S(t, α) = (1 - α)*𝒻F(t)/4

#instantaneous temperature [K]
𝒻T(t, fCO2, α=α₀) = (𝒻S(t,α) - OLR₀ + b*log(fCO2/fCO2₀))/a + T₀

#fraction of carbon in the atmosphere [-], see Mills et al., 2011
𝒻ϕ(C, h=𝐡) = 0.78*C/(C + h)

#partial pressure of CO2 [Pa]
𝒻pCO2(C, h=𝐡) = 𝒻ϕ(C,h)*C*𝛍*𝐠/𝐒ₑ

#molar concentration of CO2 [ppmv]
𝒻fCO2(C, P=𝐏₀, h=𝐡) = 1e6*𝒻pCO2(C,h)/(𝒻pCO2(C,h) + P)

#------------------------------------------------------------------------------
# integration/modeling

export integrate
function integrate(V::PowerLawDistribution, nstep::Int=10_000)
    t = LinRange(0, 1, nstep+1)
    dt = (maximum(t) - minimum(t))/nstep
    C = zeros(nstep+1)
    C[1] = 3.193e18
    for i ∈ 2:nstep+1
        C[i] = C[i-1] + sqrt(dt)*1e12*(V() - V.μ)*1e4
    end
    return t, C, 𝒻fCO2.(C)
end

export integrations
function integrations(nstep::Int=10_000, nrun::Int=10_000)::Vector{Float64}
    V = PowerLawDistribution(7.5, 250, -1.8)
    fCO2 = zeros(nrun)
    @threads for i ∈ 1:nrun
        fCO2[i] = integrate(V, nstep)[3][end]
    end
    return fCO2
end

export martingale
function martingale(; nstep::Int=100, t=1.0)
    t = LinRange(0, t, nstep+1)
    dt = (maximum(t) - minimum(t))/nstep
    x = zeros(nstep+1)
    x[1] = 0
    for i ∈ 2:nstep+1
        x[i] = x[i-1] + randn()*sqrt(dt)
    end
    return t, x
end

end