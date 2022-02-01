using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using Roots: find_zero

## calibration conditions

#target weathering rate to balance outgassing [teramole/yr]
W₀ = Vᵣ
#ocean atmosphere carbon [teramole]
C = Cᵣ
#age of solar system [Gyr]
t = 𝐭
#super strict tolerance
tol = big(1e-50)

## WHAK

k = W₀/𝒻whak(C, t, k=1.0)
println("calibrated k = $k")

## MAC

Λ = Float64(
    find_zero(
        x -> 𝒻mac(C, t, Λ=x) - W₀,
        big(1e-5),
        atol=tol,
        rtol=tol
    )
)
println("calibrated Λ = $Λ")
