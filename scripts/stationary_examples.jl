using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

𝒻fₑ(t, χ) = t |> (𝒻fCO2 ∘ χ)

step(fᵢ, fₑ, Δt, τ, σ) = fᵢ + Δt*(fₑ - fᵢ)/τ + √(Δt)*σ*randn()

#ornstein uhlenbeck process
function oup(t₁::Real, t₂::Real, nstep::Int, τ::Real, σ::Real, tspin::Real=1.0)
    @assert t₂ > t₁
    #setup
    χ = Χ() #equilibrium CO₂ over time
    t = LinRange(t₁, t₂, nstep+1)
    Δt = (t₂ - t₁)/nstep
    #spinup
    fₑ = 𝒻fₑ(t₁, χ)
    fₛ = fₑ
    tₛ = 0.0
    while tₛ < tspin
        fₛ = step(fₛ, fₑ, Δt, τ, σ)
        tₛ += Δt
    end
    #integrate
    f = zeros(Float64, nstep+1)
    f[1] = fₛ
    for i ∈ 1:nstep
        f[i] = f[i] <= 0 ? one(f[i]) : f[i]
        f[i+1] = step(f[i], 𝒻fₑ(t[i], χ), Δt, τ, σ)
    end
    println(minimum(f))
    return t, 𝒻T.(f, t)
end

##

t₁ = 4
t₂ = 4.5
nstep = 1_000_000
τ = 0.003
σ = 4000

t, T = oup(t₁, t₂, nstep, τ, σ)
n = length(t)
t, T = t[n÷2:end], T[n÷2:end]
mya = 1e3*𝒻gya.(t)
S = map(x -> (x < 280) ? x : NaN, T)

fig, ax = plt.subplots(1, 1, figsize=(4,2), constrained_layout=true)

ax.plot(mya, T, color="C1")
ax.plot(mya, S, color="C3")
ax.invert_xaxis()
ax.grid(false)
ax.plot(
    [mya[1],mya[end]],
    [280,280],
    linewidth=2.5,
    color="k",
    alpha=0.5,
    zorder=1
)
ax.set_xlim(mya[1], mya[end])
ax.set_xlabel("Time [Mya]")
ax.set_ylabel("Temperature [K]")
ax.annotate(
    raw"$T_{snow}$",
    (mya[1]-2,279.5),
    va="top",
    ha="left"
)
#fig.savefig(plotsdir("stationary_example"), dpi=500)