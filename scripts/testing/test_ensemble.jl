using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot
using Statistics
using Base.Threads: @threads

pygui(true)

##

nstep = 10_000
N = 7*250
L = 100
T = zeros(L, N)
t = zeros(L, N)
t₁ = 2.5
idx = Int64.(round.(range(1, nstep, L)))
for i in 1:N
    tᵢ, Cᵢ, Vᵢ = simulate(
        initparams(
            μ=Vᵣ,
            τ=1e7,
            σ=1e-4,
            Vₘ=0,
        ),
        t₁=t₁,
        C₁=𝒻Cₑ(t₁),
        𝒻W=(C,t)->𝒻whak(C, t, β=0),
        nstep=nstep
    )
    Tᵢ = @. 𝒻T(𝒻fCO2(Cᵢ), tᵢ)
    T[:,i] = Tᵢ[idx]
    t[:,i] = tᵢ[idx]
end

##

fig, ax = plt.subplots(1,1)
for q ∈ (0.001, 0.01, 0.1, 0.25)
    ax.plot(
        t[:,1],
        [quantile(T[i,:], 1 - q) for i in 1:L],
        color="k",
        alpha=0.2
    )
    ax.plot(
        t[:,1],
        [quantile(T[i,:], q) for i in 1:L],
        color="k",
        alpha=0.2
    )
end
ax.plot(
    t[:,1],
    [median(T[i,:]) for i in 1:L],
    color="k"
)
ax.grid(false)
