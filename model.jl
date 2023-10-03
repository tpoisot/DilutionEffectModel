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
include(joinpath("lib", "helpers.jl"))

# Generate a dataframe for simulations
bank = DataFrame()
centers = collect(0:10:100)
barycenters = [(c1, c2, c3) for c1 in centers for c2 in centers for c3 in centers]
filter!(b -> isequal(100)(sum(b)), barycenters)

for b in barycenters
    push!(bank, (
        id=uuid4(),
        mutualism=b[1],
        competition=b[2],
        predation=b[3],
    ))
end

progressbar = Progress(size(bank, 1));

results = [DataFrame() for thr in 1:Threads.nthreads()]

Threads.@threads for parameters in eachrow(bank)

    sim_id = parameters.id
    intprop = [parameters.mutualism, parameters.competition, parameters.predation]

    for replicate in 1:150
        Sᵢ, Iᵢ, Hᵢ = onesim(S, intprop)
        if ~isempty(Hᵢ)
            repl_id = uuid4()
            push!(results[Threads.threadid()], (
                parameters=sim_id,
                replicate=repl_id,
                S=sum(Sᵢ),
                I=sum(Iᵢ),
                H=sum(Hᵢ),
                prevalence=sum(Iᵢ) / sum(Hᵢ),
                diversity=shannon(Hᵢ),
                richness=length(Hᵢ)
            ))
        end
    end

    next!(progressbar)
end

results = vcat(results...)
data = dropmissing!(leftjoin(results, bank, on=:parameters => :id))

final = combine(
    groupby(data, :parameters),
    :prevalence => mean => :prevalence,
    :prevalence => safe_std => :prevalence_std,
    :diversity => mean => :diversity,
    :richness => median => :richness,
    [:diversity, :prevalence] => safe_cor => :correlation,
    :mutualism => first => :mutualism,
    :competition => first => :competition,
    :predation => first => :predation,
    nrow => :n
)

begin
    fig = Figure(; resolution=(800, 800))
    ax = Axis(fig[1, 1])

    mut = final.mutualism / 100
    cmp = final.competition / 100
    prd = final.predation / 100

    ternaryaxis!(
        ax;
        labelx="Mutualism",
        labely="Competition",
        labelz="Predation"
    )

    divpal = ColorSchemes.diverging_bwg_20_95_c41_n256
    linpal = ColorSchemes.linear_gow_65_90_c35_n256

    zval = final.correlation
    cpal = minimum(zval) < 0.0 ? divpal : linpal
    mval = minimum(zval) < 0.0 ? (-1., 1.) : extrema(zval)
    point_color = get(cpal, zval, mval)

    tplot = ternaryscatter!(
        ax,
        mut, cmp, prd,
        color=point_color,
        marker=:hexagon,
        markersize=35,
    )

    Colorbar(fig[end+1, 1], colormap=cpal, limits=mval, vertical=false)

    xlims!(ax, -0.2, 1.2)
    ylims!(ax, -0.3, 1.1)
    hidedecorations!(ax)
    fig
end
