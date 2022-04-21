using DrWatson
@quickactivate "Random Volcanic Climate"
push!(LOAD_PATH, srcdir())
using RandomVolcanicClimate
using Distributions
using PyPlot

pygui(true)

##

#sets axis y limits to the same values
function sharey(axs)
    ymin = minimum(ax->ax.get_ylim()[1], axs)
    ymax = maximum(ax->ax.get_ylim()[2], axs)
    foreach(ax->ax.set_ylim(ymin, ymax), axs)
    nothing
end

##

cmap = plt.cm.get_cmap("cool")

##

#time slice selections [Gyr]
t = [4.1, 4.3, 4.5]
#fCO2 standard deviations [ppm]
σ = [50, 100, 150, 200]

#initialize axes
fig, axs = plt.subplots(2, length(t), figsize=(7,5))
#convert times to Mya
mya = @. Int(round(1e3*(4.5 - t)))

for i in 1:size(axs,2)
    Cₑ = 𝒻Cₑ(t[i]) #equilibrium carbon
    fₑ = 𝒻fCO2(Cₑ) #equilibrium CO2
    for j ∈ 1:length(σ)
        #construct the normal distribution
        N = Normal(fₑ, σ[j])
        #create a range of fCO2 values between prescribed quantiles
        f = LinRange(
            max(quantile(N, 1e-6), 0),
            quantile(N, 1 - 1e-6),
            10_000
        )
        #plot the co2 pdf
        axs[1,i].plot(
            pdf.(N, f),
            f,
            color=cmap((j-1)/(length(σ)-1)),
            linewidth=1.75,
            label="σ=$(σ[j])"
        )
        #compute temperature range corresponding to fCO2 range
        T = 𝒻T.(f, t[i])
        #plot the same pdf against temperature coordinates
        axs[2,i].plot(
            pdf.(N, f),
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
    axs[1,j].set_xticks([0,axs[1,j].get_xticks()[end]])
    axs[1,j].set_xticklabels([])
    axs[2,j].set_xticks([0,axs[2,j].get_xticks()[end]])
    axs[2,j].set_ylim(275, axs[2,j].get_ylim()[2])
    axs[1,j].set_title("$(mya[j]) Mya")
end
#foreach(axs) do ax
#    ax.set_xlim(-0.01,1.04)
#end
fig.supxlabel("Probability Density")
plt.tight_layout()
fig.savefig(plotsdir("T_pdfs"), dpi=500)

## probability densities of snowball temperatures through time

function Psnow(σ, t, Tsnow)
    Cₑ = 𝒻Cₑ(t) #Tmole
    fₑ = 𝒻fCO2(Cₑ) #ppm
    N = Normal(fₑ, σ)
    Csnow = 𝒻Cₑ(t, Tsnow)
    fsnow = 𝒻fCO2(Csnow)
    pdf(N, fsnow)
end

##

fig, ax = plt.subplots(1, 1, figsize=(4,3.5))
t = LinRange(3.5, 4.5, 1001)
σ = [50, 100, 150, 200]
Tsnow = 280
for i ∈ 1:length(σ)
    ax.semilogy(
        gya.(t),
        Psnow.(σ[i], t, Tsnow),
        label="σ=$(σ[i])",
        color=cmap((i-1)/(length(σ)-1)),
        linewidth=2
    )
end
ax.invert_xaxis()
ax.set_ylim(1e-16, 1)
ax.legend()
ax.set_xlabel("Time [Gya]")
ax.set_title("Snowball (T = $Tsnow K) Probability Density")
fig.tight_layout()
fig.savefig(plotsdir("Tsnow_probabiliy_densities"), dpi=500)
