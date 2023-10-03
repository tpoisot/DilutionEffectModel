using CairoMakie
using DifferentialEquations
using TernaryDiagrams
using Distributions
using StatsBase
using LinearAlgebra
import ColorSchemes
using ProgressMeter
using Statistics
using UUIDs
using DataFrames
import CSV

# Load some functions
include(joinpath("lib", "parameters.jl"))
include(joinpath("lib", "interactions.jl"))
include(joinpath("lib", "ode.jl"))
include(joinpath("lib", "onesim.jl"))

# Generate a dataframe for simulations
simulation_bank = DataFrame()
centers = collect(0:4:100)
barycenters = [(c1, c2, c3) for c1 in centers for c2 in centers for c3 in centers]
filter!(b -> isequal(100)(sum(b)), barycenters)

for b in barycenters
    push!(simulation_bank, (
        id = uuid4(),
        mutualism = b[1],
        competition = b[2],
        predation = b[3],
    ))
end

progressbar = Progress(size(simulation_bank, 1));

simulations_results = [DataFrame() for thr in 1:Threads.nthreads()]

Threads.@threads for i in 1:size(simulation_bank, 1)
    
    sim_id = simulation_bank.id[i]
    intprop = range(simulation_bank.mutualism[i], simulation_bank.competition[i], simulation_bank.predation[i])
    
    for replicate in 1:5
        Sᵢ, Iᵢ, Hᵢ = onesim(S, intprop)
        if ~isempty(Hᵢ)
            repl_id = uuid4()
            push!(simulations_results[Threads.threadid()], (
                parameters = sim_id,
                replicate = repl_id,
                S = sum(Sᵢ),
                I = sum(Iᵢ),
                H = sum(Hᵢ),
                prevalence = sum(Iᵢ)/sum(Hᵢ),
                diversity = shannon(Hᵢ),
                richness = length(Hᵢ)
            ))
        end
    end
    
    next!(progressbar)
end

begin
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

    tplot = ternarycontour!(
        ax,
        cmp, mut, prd,
        richness,
        colormap = :RdYlGn, levels=10
    )

    Colorbar(fig[1,end+1], tplot)

    xlims!(ax, -0.2, 1.2)
    ylims!(ax, -0.3, 1.1)
    hidedecorations!(ax)
    fig
end