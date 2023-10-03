using CairoMakie
using DifferentialEquations
using TernaryDiagrams
using Distributions
using StatsBase
using LinearAlgebra
import ColorSchemes
using ProgressMeter

# Load some functions
include(joinpath("lib", "parameters.jl"))
include(joinpath("lib", "interactions.jl"))
include(joinpath("lib", "ode.jl"))

function onepopislow(u, t, integrator)
    return length(u) - count(u .< 10eps())
end
function extinguishpop!(integrator)
    for i in eachindex(integrator.u)
        if integrator.u[i] < 10eps()
            integrator.u[i] = 0.0
        end
    end
end
cb = ContinuousCallback(onepopislow, extinguishpop!)

function shannon(x)
    p = x ./ sum(x)
    return - sum(p .* log.(p))
end

range(x...) = x ./ sum(x)

# Centers for the populations
centers = collect(10:5:90)
barycenters = [(c1, c2, c3) for c1 in centers for c2 in centers for c3 in centers]
filter!(b -> isequal(110)(sum(b)), barycenters)

progressbar = Progress(length(barycenters));

diversity = zeros(Float64, length(barycenters))
prevalences = zeros(Float64, length(barycenters))
richness = zeros(Float64, length(barycenters))

Threads.@threads for i in eachindex(barycenters)
    
    intprop = range(barycenters[i]...) # Mutualism, Competition, Predation
    
    interactions!(A, Co, K, σₓ, σᵢ, intprop)
    
    β = transmissionmatrix(A, δ)
    
    M = A[2,:,:]
    C = A[3,:,:]
    P = A[4,:,:]
    
    u₀ = 10.0rand(2S)
    p = (S, r, δ, h, β, M, P, C, ρ, ν)

    prob = ODEProblem(densitydependent, u₀, (0.0, timesteps), p)
    sol = solve(prob, callback=cb)

    SI = sol[end]
    Sₜ = SI[1:S]
    Iₜ = SI[(S+1):end]
    Hₜ = Sₜ .+ Iₜ
    
    remain = findall((Sₜ .> 10eps()).*(Iₜ .> 10eps()))
    Sₜ = Sₜ[remain]
    Iₜ = Iₜ[remain]
    Hₜ = Sₜ .+ Iₜ

    diversity[i] = shannon(Hₜ)
    richness[i] = length(remain)
    prevalences[i] = sum(Iₜ) / sum(Hₜ)

    next!(progressbar)
end

fig = Figure();
ax = Axis(fig[1, 1]);

x1 = map(x -> x[1]/sum(x), barycenters)
x2 = map(x -> x[2]/sum(x), barycenters)
x3 = map(x -> x[3]/sum(x), barycenters)

cs = [get(ColorSchemes.lapaz, w, extrema(richness)) for w in richness]

ternaryaxis!(ax);
ternaryscatter!(
    ax,
    x1, x2, x3,
    color = cs,
    marker = :circle,
    markersize = 20,
)

xlims!(ax, -0.2, 1.2)
ylims!(ax, -0.3, 1.1)
hidedecorations!(ax)
fig