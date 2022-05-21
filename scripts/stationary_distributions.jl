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

#computes the probability of the fCO₂ value where snowball temperature is achieved
function Psnow(σ, t, Tsnow)
    Cₑ = 𝒻Cₑ(t) #Tmole
    fₑ = 𝒻fCO2(Cₑ) #ppm
    N = Normal(fₑ, σ) #fCO2 distribution
    Csnow = 𝒻Cₑ(t, Tsnow)
    fsnow = 𝒻fCO2(Csnow)
    pdf(N, fsnow)
end

#computes the temperature at a specified fCO₂ quantile
function Tquantile(σ, t, q)
    Cₑ = 𝒻Cₑ(t) #Tmole
    fₑ = 𝒻fCO2(Cₑ) #ppm
    N = Normal(fₑ, σ) #fCO2 distribution
    f = quantile(N, q) #fCO2 value of specified quantile
    𝒻T(max(f,0.0), t) #temperature at the fCO2 quantile
end

##

#time slice selections [Gyr]
t = [4, 4.25, 4.5]
#fCO2 standard deviations [ppm]
σ = [50, 100, 150, 200]

#initialize axes
fig, axs = plt.subplots(2, length(t), figsize=(7,4.5), constrained_layout=true)
#colormap selection
cmap = plt.cm.get_cmap("cool")
#convert times to Mya
mya = @. Int(round(1e3*𝒻gya(t)))

for i in 1:size(axs,2)
    Cₑ = 𝒻Cₑ(t[i]) #equilibrium carbon
    fₑ = 𝒻fCO2(Cₑ) #equilibrium CO2
    for j ∈ 1:length(σ)
        #construct the normal distribution
        N = Normal(fₑ, σ[j])
        #create a range of fCO2 values between prescribed quantiles
        f = LinRange(
            max(quantile(N, 1e-8), 0),
            quantile(N, 1 - 1e-8),
            1000
        )
        #plot the co2 pdf
        axs[1,i].plot(
            pdf.(N, f)/pdf(N, fₑ),
            f,
            color=cmap((j-1)/(length(σ)-1)),
            linewidth=1.75,
            label="σ=$(σ[j])"
        )
        #compute temperature range corresponding to fCO2 range
        T = 𝒻T.(f, t[i])
        #plot the same pdf against temperature coordinates
        axs[2,i].plot(
            pdf.(N, f)/pdf(N, fₑ),
            T,
            color=cmap((j-1)/(length(σ)-1)),
            linewidth=1.75,
            zorder=0,
            label="σ=$(σ[j])"
        )
    end
    axs[1,i].plot([0,1], [fₑ,fₑ], "k", linewidth=0.75, alpha=0.75)
    axs[2,i].plot([0,1], [288,288], "k", linewidth=0.75, alpha=0.75)
    axs[1,i].annotate(
        raw"$\mu_{fCO_2}=$" * (fₑ |> round |> Int |> string),
        (1,fₑ+250),
        va="bottom",
        ha="right"
    )
    fig.add_artist(
        matplotlib.patches.ConnectionPatch(
            xyA=(0.5,-0.05),
            xyB=(0.5,0.9),
            coordsA="axes fraction",
            coordsB="axes fraction",
            axesA=axs[1,i],
            axesB=axs[2,i],
            arrowstyle="->",
            lw=2
        )
    )
end
axs[2,1].legend()
axs[1,1].set_ylabel("CO₂ [ppm]")
axs[2,1].set_ylabel("Temperature [K]")
sharey(axs[1,:])
sharey(axs[2,:])
foreach(axs) do ax
    ax.grid(false)
    ax.set_xlim(0, 1.02)
    ax.set_xticks([0,1])
end
for j in 2:size(axs,2)
    axs[1,j].set_yticklabels([])
    axs[2,j].set_yticklabels([])
end
for j in 1:size(axs,2)
    axs[1,j].set_xticklabels([])
    axs[1,j].set_ylim(0, axs[1,j].get_ylim()[2])
    axs[2,j].set_ylim(275, axs[2,j].get_ylim()[2])
    axs[2,j].plot(
        axs[2,j].get_xlim(),
        [280,280],
        linewidth=2.5,
        color="k",
        alpha=0.5,
        zorder=-1
    )
    axs[1,j].set_title("$(mya[j]) Mya")
end
axs[2,end].annotate(
    raw"$T_{snow}$",
    (axs[2,end].get_xlim()[2],280.5),
    va="bottom",
    ha="right"
)
fig.supxlabel("Normalized Probability Density")
fig.savefig(plotsdir("stationary_distributions"), dpi=500)

## probability densities of snowball temperatures through time

fig, axs = plt.subplots(1, 2, figsize=(7,3.5), constrained_layout=true)
t = LinRange(4, 4.5, 1001)
gya = 𝒻gya.(t)*1e3

σ = [50, 100, 150, 200]
Tsnow = 280.0
cmap = plt.cm.get_cmap("cool")
for i ∈ 1:length(σ)
    color = cmap((i-1)/(length(σ)-1))
    y = Psnow.(σ[i], t, Tsnow)
    axs[1].semilogy(gya, y, color=color, linewidth=2)
    y = Tquantile.(σ[i], t, 1e-2)
    axs[2].plot(gya, y, label="σ=$(σ[i])", color=color, linewidth=2)
end

axs[1].set_xlim(minimum(gya), maximum(gya))
axs[1].set_ylim(1e-16, 1)
axs[1].invert_xaxis()
axs[1].set_xlabel("Time [Mya]")
axs[1].set_title("Snowball Probability Density")

axs[2].set_xlim(minimum(gya), maximum(gya))
axs[2].set_ylim(275, axs[2].get_ylim()[2])
axs[2].invert_xaxis()
axs[2].legend()
axs[2].set_xlabel("Time [Mya]")
axs[2].set_title("Temperature [K] of 1ˢᵗ CO₂ Percentile")

fig.savefig(plotsdir("stationary_snowball_probabilities"), dpi=500)
