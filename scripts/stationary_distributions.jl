using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using Distributions
using PyPlot

pygui(true)

##

function sharey(axs)
    ymin = minimum(ax->ax.get_ylim()[1], axs)
    ymax = maximum(ax->ax.get_ylim()[2], axs)
    foreach(ax->ax.set_ylim(ymin, ymax), axs)
    nothing
end

##

cmap = plt.cm.get_cmap("cool")

##

t = [4.2, 4.3, 4.4, 4.5]
fig, axs = plt.subplots(2, length(t), figsize=(7,5))
mya = @. Int(round(1e3*(4.5 - t)))
σ = [50, 100, 150, 200]

for i in 1:size(axs,2)
    Cₑ = 𝒻Cₑ(t[i]) #equilibrium carbon
    fₑ = 𝒻fCO2(Cₑ) #equilibrium CO2
    for j ∈ 1:length(σ)
        #construct the normal distribution
        N = Normal(fₑ, σ[j])
        #determine a co2 range
        f = LinRange(
            max(quantile(N, 1e-4), 0),
            quantile(N, 1 - 1e-4),
            10_000
        )
        #plot the co2 pdf
        axs[1,i].plot(
            pdf.(N, f)/pdf(N, fₑ),
            f,
            color=cmap((j-1)/(length(σ)-1)),
            linewidth=1.75,
            label="σ=$(σ[j])"
        )
        #compute resulting temperature range
        T = 𝒻T.(f, t[i])
        #plot T pdf
        axs[2,i].plot(
            pdf.(N, f)/pdf(N, fₑ),
            T,
            color=cmap((j-1)/(length(σ)-1)),
            linewidth=1.75,
        )
    end
    
end
axs[1,end].legend()
axs[1,1].set_ylabel("CO₂ [ppm]")
axs[2,1].set_ylabel("Temperature [K]")
sharey(axs[1,:])
sharey(axs[2,:])
for j in 2:size(axs,2)
    axs[1,j].set_yticklabels([])
    axs[2,j].set_yticklabels([])
end
for j in 1:size(axs,2)
    axs[1,j].set_xticks([0,1])
    axs[1,j].set_xticklabels([])
    axs[2,j].set_xticks([0,1])
    axs[2,j].set_ylim(275, axs[2,j].get_ylim()[2])
    axs[1,j].set_title("$(mya[j]) Mya")
end
foreach(axs) do ax
    ax.set_xlim(0,1.04)
end
fig.supxlabel("Normalized Probability Distribution")
plt.tight_layout()
fig.savefig(plotsdir("T_pdfs"), dpi=500)
