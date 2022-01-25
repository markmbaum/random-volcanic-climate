export PowerLawDistribution, ∫

struct PowerLawDistribution
    x₁::Float64
    x₂::Float64
    C::Float64
    γ::Float64
    μ::Float64
    σ²::Float64
end

function Base.show(io::IO, D::PowerLawDistribution)
    @unpack x₁, x₂, γ, μ, v = D
    println(io, "PowerLawDistribution")
    println(io, "  x₁ = $x₁")
    println(io, "  x₂ = $x₂")
    println(io, "  γ = $γ")
    println(io, "  μ = $μ")
    print(io, "  v = $v")
end

function powerlawmean(x₁, x₂, γ)
    x₁ == x₂ && return x₁
    g = γ + 1
    h = γ + 2
    g*(x₂^h - x₁^h)/(h*(x₂^g - x₁^g))
end

function PowerLawDistribution(μ::Real, x₂::Real, γ::Real)
    #@assert γ < -1 "gamma must be < -1"
    @assert !(γ == Int(round(γ))) "gamma can't be a round integer"
    #use extended arithmetic to get an exact value
    x₀ = big(μ/2)
    n::Int64 = 0
    while powerlawmean(x₀, x₂, γ) > μ
        x₀ /= 10
        n += 1
        n > 100 && error("failure to construct power law distribution")
    end
    tol = big(1e-30)
    x₁ = Float64(
        find_zero(
            x -> powerlawmean(x, big(x₂), big(γ)) - big(μ),
            (x₀, big(x₂)),
            atol=tol,
            rtol=tol,
            xatol=tol,
            xrtol=tol
        )
    )
    #normalization coefficient
    g = γ + 1
    C = g/(x₂^g - x₁^g)
    #variance
    h = γ + 2
    k = γ + 3
    σ² = (g/(x₂^g - x₁^g))*((x₂^k - x₁^k)/k - 2μ*(x₂^h - x₁^h)/h + μ^2*(x₂^g - x₁^g)/g)
    #construct
    PowerLawDistribution(x₁, x₂, C, γ, μ, σ²)
end

function sample(D::PowerLawDistribution, U::Float64)::Float64
    @unpack x₁, x₂, γ = D
    g = γ + 1
    ((x₂^g - x₁^g)*U + x₁^g)^(1/g)
end

(D::PowerLawDistribution)()::Float64 = sample(D, rand())

(D::PowerLawDistribution)(rng::AbstractRNG)::Float64 = sample(D, rand(rng))

(D::PowerLawDistribution)(n::Int) = map(u->sample(D,u), rand(n))

(D::PowerLawDistribution)(rng::AbstractRNG, n::Int) = map(u->sample(D,u), rand(rng, n))

function ∫(D::PowerLawDistribution)
    @unpack x₁, x₂, γ, C = D
    hquadrature(x -> C*x^γ, x₁, x₂)[1]
end