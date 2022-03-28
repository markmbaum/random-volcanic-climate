using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using Base.Threads: @threads, nthreads
using PyPlot

pygui(true)

##

𝒻W(C,t) = 𝒻whak(C, t, β=0)
t₁ = 2.5

##

N = 500*nthreads()
tsnow = zeros(N)

@threads for i in 1:N
    t, C, V = simulate(
        initparams(
            μ=Vᵣ,
            τ=1e6,
            σ=4e-3,
            Vₘ=0
        ),
        t₁=t₁,
        C₁=𝒻Cₑ(t₁),
        𝒻W=𝒻W,
        nstep=10_000
    )
    T = C2T.(C, t)
    tsnow[i] = snowballtime(t, T)
end

hist(tsnow);