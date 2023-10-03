using CairoMakie
using DifferentialEquations
using TernaryDiagrams
using Distributions
using StatsBase
using LinearAlgebra
import ColorSchemes
using ProgressMeter
using Statistics

# Load some functions
include(joinpath("lib", "parameters.jl"))
include(joinpath("lib", "interactions.jl"))
include(joinpath("lib", "ode.jl"))
include(joinpath("lib", "onesim.jl"))

function shannon(x)
    p = x ./ sum(x)
    return - sum(p .* log.(p))
end

range(x...) = x ./ sum(x)

# Centers for the populations
centers = collect(0:10:100)
barycenters = [(c1, c2, c3) for c1 in centers for c2 in centers for c3 in centers]
filter!(b -> isequal(100)(sum(b)), barycenters)

progressbar = Progress(length(barycenters));

diversity = zeros(Float64, length(barycenters))
prevalences = zeros(Float64, length(barycenters))
richness = zeros(Float64, length(barycenters))
correlations = zeros(Float64, length(barycenters))

Threads.@threads for i in eachindex(barycenters)
    
    intprop = range(barycenters[i]...) # Mutualism, Competition, Predation
    
    #Sₜ, Iₜ, Hₜ = onesim(S, intprop)
    outputs = [onesim(S, intprop) for _ in 1:20]
    filter!(x -> !isempty(x[3]), outputs)

    div = [shannon(output[3]) for output in outputs]
    ric = [length(output[3]) for output in outputs]
    prv = [sum(output[2])/sum(output[3]) for output in outputs]
    cdp = isempty(outputs) ? 0.0 : cor(div, prv)

    diversity[i] = mean(div)
    richness[i] = mean(ric)
    prevalences[i] = mean(prv)
    correlations[i] = isnan(cdp) ? 0.0 : cdp

    next!(progressbar)
end

fig = Figure();
ax = Axis(fig[1, 1]);

mut = map(x -> x[1]/sum(x), barycenters)
cmp = map(x -> x[2]/sum(x), barycenters)
prd = map(x -> x[3]/sum(x), barycenters)

ternaryaxis!(
    ax;
    labelx = "Mutualism",
    labely = "Competition",
    labelz = "Predation",
);

tplot = ternarycontourf!(
    ax,
    cmp, mut, prd,
    prevalences,
    colormap = :dense, levels=10
)

Colorbar(fig[1,end+1], tplot)

xlims!(ax, -0.2, 1.2)
ylims!(ax, -0.3, 1.1)
hidedecorations!(ax)
fig