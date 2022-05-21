using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate

##

for ensname ∈ ("wide", "deep")
    fn = datadir("sims", "ensemble_"*ensname*".jld2")
    V = (fn |> loadensemble)[:R][:V]
    Tₑ = 𝒻Tₑ(𝒻whak, V)
    fn = datadir("exp_pro", "Te_"*ensname*".jld2")
    safesave(fn, Dict("Tₑ"=>Tₑ))
end
