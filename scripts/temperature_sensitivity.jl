using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using ForwardDiff: derivative
using RandomVolcanicClimate
using PyPlot

pygui(true)

##

𝒻fₑ(t::Real, χ) = t |> (𝒻fCO2 ∘ χ)

𝒻fₑ(t::AbstractVector, χ) = map(x -> 𝒻fₑ(x, χ), t)

∂T(f, t=𝐭) = derivative(x -> 𝒻T(x, t), f)

##

t = LinRange(2.5, 4.5, 101)
f = 𝒻fₑ(t, Χ())
∂ = ∂T.(f, t)
gya = 𝒻gya.(t)

##

fig, axf = plt.subplots(1, 1, figsize=(4,2.4), constrained_layout=true)
axf.set_xlabel("Time [Gya]")
axf.tick_params(axis="y", which="both", labelcolor="C1", colors="C1")
axf.set_ylabel("χ [ppm]", color="C1")
axf.semilogy(gya, f, color="C1")
axf.grid(false)
axf.invert_xaxis()

ax∂ = axf.twinx()
ax∂.tick_params(axis="y", labelcolor="k")
ax∂.set_ylabel("∂T/∂f [K/ppm]", rotation=270, va="bottom")
ax∂.plot(gya, ∂, color="k")
ax∂.grid(false)
ax∂.set_ylim(0, ax∂.get_ylim()[2])
ax∂.spines["right"].set_visible(true)
ax∂.spines["top"].set_visible(true)

axf.spines["left"].set_color("C1")
ax∂.spines["left"].set_color("C1")

axf.set_title("Temperature Sensitivity to CO₂ Perturbation")
fig.savefig(plotsdir("temperature_sensitivity"), dpi=500)
plt.close(fig)
