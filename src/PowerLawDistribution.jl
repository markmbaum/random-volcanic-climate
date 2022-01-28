module PowerLawDistribution

using Distributions
using UnPack
using Roots
using Random: AbstractRNG
using Cubature: hquadrature

import Base: rand
import Statistics: mean, var

export PowerLaw, ∫, mean, var, rand

struct PowerLaw <: Sampleable{Univariate,Continuous}
    x₁::Float64
    x₂::Float64
    C::Float64
    γ::Float64
    μ::Float64
    σ²::Float64
end

function Base.show(io::IO, D::PowerLaw)
    @unpack x₁, x₂, γ, μ, σ² = D
    println(io, "PowerLaw")
    println(io, "  x₁ = $x₁")
    println(io, "  x₂ = $x₂")
    println(io, "  γ = $γ")
    println(io, "  μ = $μ")
    print(io, "  σ² = $σ²")
end

function plmean(x₁, x₂, γ)
    x₁ == x₂ && return x₁
    g = γ + 1
    h = γ + 2
    μ = g*(x₂^h - x₁^h)/(h*(x₂^g - x₁^g))
    return μ
end

function plvar(x₁, x₂, γ)
    x₁ == x₂ && return zero(x₁)
    g = γ + 1
    h = γ + 2
    k = γ + 3
    μ = plmean(x₁, x₂, γ)
    σ² = (g/(x₂^g - x₁^g))*((x₂^k - x₁^k)/k - 2μ*(x₂^h - x₁^h)/h + μ^2*(x₂^g - x₁^g)/g)
    return σ²
end

function plnorm(x₁, x₂, γ)
    g = γ + 1
    C = g/(x₂^g - x₁^g)
    return C
end

function plx₁(μ, x₂, γ)
    x₀ = float(μ)
    n::Int64 = 0
    while plmean(x₀, x₂, γ) > μ
        x₀ /= 10
        n += 1
        n > 100 && error("x₀=$x₀, failure to construct power law distribution")
    end
    tol = 1e-12
    find_zero(
        x -> plmean(x, x₂, γ) - μ,
        (x₀, x₂),
        atol=tol,
        rtol=tol,
        xatol=tol,
        xrtol=tol
    )
end

function ∫(D::PowerLaw)
    @unpack x₁, x₂, γ, C = D
    hquadrature(x -> C*x^γ, x₁, x₂)[1]
end

function PowerLaw(μ::Real, x₂::Real, γ::Real)
    @assert γ < -1 "gamma must be < -1"
    @assert !(γ ∈ (-1,-2,-3)) "gamma can't be -1, -2, or -3"
    #find the lower boundary of the distribution
    x₁ = (μ == x₂) ? μ : plx₁(μ, x₂, γ)
    #normalization coefficient
    C = plnorm(x₁, x₂, γ)
    #variance
    σ² = plvar(x₁, x₂, γ)
    #construct
    PowerLaw(x₁, x₂, C, γ, μ, σ²)
end

#--------------------------------------
# sampling

mean(D::PowerLaw) = D.μ

var(D::PowerLaw) = D.σ²

#--------------------------------------
# sampling

function sample(D::PowerLaw, 𝒰)
    @unpack x₁, x₂, γ = D
    g = γ + 1
    ((x₂^g - x₁^g)*𝒰 + x₁^g)^(1/g)
end

rand(D::PowerLaw) = sample(D, Base.rand())

rand(D::PowerLaw, n::Int) = map(i->rand(D), 1:n)

rand(rng::AbstractRNG, D::PowerLaw) = sample(D, rand(rng))

rand(rng::AbstractRNG, D::PowerLaw, n::Int) = map(i->rand(rng, D), 1:n)

end
