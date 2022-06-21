using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

𝒻fₑ(t::Real, χ) = t |> (𝒻fCO2 ∘ χ)

𝒻fₑ(t::AbstractVector, χ) = map(x -> 𝒻fₑ(x, χ), t)

step(fᵢ, fₑ, Δt, τ, σ) = fᵢ + Δt*(fₑ - fᵢ)/τ + √(Δt)*σ*randn()

#ornstein uhlenbeck process
function OUP(t₁::Real, t₂::Real, nstep::Int, τ::Real, σ::Real, tspin::Real=1.0)
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
    return t, f, 𝒻T.(f, t)
end

##

function annotate(ax, s)
    ax.annotate(
        s,
        (0.01,0.99),
        xycoords="axes fraction",
        va="top",
        fontsize=12,
        fontweight="bold"
    )
    nothing
end

function Tsnow(ax, a, b)
    ax.plot(
        [a, b],
        [280,280],
        linewidth=2.5,
        color="k",
        alpha=0.5,
        zorder=1
    )
    nothing
end

##

t₁ = 4
t₂ = 4.5
nstep = 1_000_000
τ = 0.003
σ = 3000

t, f, T = OUP(t₁, t₂, nstep, τ, σ)
n = length(t)
t, f, T = t[n÷2:end], f[n÷2:end], T[n÷2:end]
mya = 1e3*𝒻gya.(t)
S = map(x -> (x < 280) ? x : NaN, T)

fig, axs = plt.subplots(
    2, 2,
    figsize=(7,3.5),
    constrained_layout=true
)


χ = Χ()
axs[1].plot(mya, f, color="C1", linewidth=0.75)
axs[1].plot(mya, 𝒻fₑ(t, χ), color="k", alpha=0.5)
axs[1].invert_xaxis()
axs[1].set_xlim(mya[1], mya[end])
axs[1].set_ylabel("CO₂ [ppm]")
axs[1].annotate(
    "χ",
    (mya[end], 𝒻fₑ(t[end], χ)),
    va="top",
    ha="right",
    fontsize=12,
    fontweight="bold"
)
annotate(axs[1], "a")


axs[2].plot(mya, T, color="C2", linewidth=0.75)
axs[2].plot(mya, S, color="dodgerblue", linewidth=0.75)
axs[2].invert_xaxis()
axs[2].plot(
    [mya[1],mya[end]],
    [280,280],
    linewidth=2.5,
    color="k",
    alpha=0.5,
    zorder=1
)
axs[2].set_xlim(mya[1], mya[end])
axs[2].set_ylabel("Temperature [K]")
Tsnow(axs[2], mya[1]-2, 279.5)
annotate(axs[2], "c")


_, idx = findmin(T)
δ = (0.0025*length(t)/(t[end] - t[1])) |> round |> Int
sl = idx-δ:idx+δ


axs[3].plot(mya[sl], f[sl], color="C1", linewidth=0.75)
axs[3].plot(mya[sl], 𝒻fₑ(t[sl], χ), color="k", alpha=0.5)
axs[3].invert_xaxis()
axs[3].set_xlim(mya[sl][1], mya[sl][end])
axs[3].annotate(
    "χ",
    (mya[sl][end], 𝒻fₑ(t[sl][end], χ)),
    va="top",
    ha="right",
    fontsize=12,
    fontweight="bold"
)
annotate(axs[3], "b")


axs[4].plot(mya[sl], T[sl], color="C2", linewidth=0.75)
axs[4].plot(mya[sl], S[sl], color="dodgerblue", linewidth=0.75)
axs[4].invert_xaxis()
axs[4].plot(
    [mya[sl][1],mya[sl][end]],
    [280,280],
    linewidth=2.5,
    color="k",
    alpha=0.5,
    zorder=1
)
axs[4].set_xlim(mya[sl][1], mya[sl][end])
Tsnow(axs[4], mya[sl][1], 279.5)
axs[4].set_ylim(axs[2].get_ylim())
annotate(axs[4], "d")
foreach(axs) do ax
    ax.grid(false)
end
fig.supxlabel("Time [Mya]")

fig.savefig(plotsdir("stationary_example"), dpi=500)
plt.close(fig)