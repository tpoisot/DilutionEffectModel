using DifferentialEquations
using Distributions
using StatsBase
using LinearAlgebra
using UUIDs
using DataFrames
import CSV
import Random
using DelimitedFiles

# Prepare to load things
input_path = joinpath(dirname(Base.active_project()), "inputs")
if ~isfile(joinpath(input_path, "sets.csv"))
    include("00_generate_simulation_parameters.jl")
end

# Load the additional functions
lib_path = joinpath(dirname(Base.active_project()), "lib")
lib_files = readdir(lib_path, join=true)
for lib_file in lib_files
    include(lib_file)
end

# Prepare the outputs folder
output_path = joinpath(dirname(Base.active_project()), "outputs")
if ~ispath(output_path)
    mkpath(output_path)
end

# Read the parameter sets
parameter_sets = DataFrame(CSV.File(joinpath(input_path, "sets.csv")))

# If the job is not running from SLURM using the full parameter set, this will sample one parameter at random
parameter_set_id = get(ENV, "SLURM_ARRAY_TASK_ID", rand(1:size(parameter_sets, 1)))
if !(typeof(parameter_set_id) <: Number)
    # The environmental variables are strings, so we need to convert
    parameter_set_id = parse(Int, parameter_set_id)
end
parameters = parameter_sets[parameter_set_id, :]

# Seed the RNG for the node
rng = Random.seed!(parameters.seed)

# Package the parameters for the run
model_parameters = (
    timesteps=parameters.time,
    Co=parameters.connectance,
    K=parameters.carrying_capacity,
    σₓ=parameters.sigma_interaction,
    σᵢ=parameters.sigma_infection,
    δ=parameters.regulation,
    r=fill(1.0, parameters.richness),
    ν=fill(parameters.virulence, parameters.richness),
    ρ=fill(parameters.recovery, parameters.richness),
    h=parameters.half_saturation,
    density=parameters.density,
)

# Interaction proportions
interaction_parameters = (
    parameters.mutualism, parameters.competition, parameters.predation
)

# Generate the UUID4s for the replicates
replicates_id = [uuid4(rng) for i in 1:parameters.replicates]

# Main loop using multi-threading
Threads.@threads for replicate in 1:parameters.replicates

    # Each replicate has a unique id
    this_replicate_id = replicates_id[replicate]

    # The simulation will return a time series, and a series of matrices
    output = onesim(parameters.richness, interaction_parameters, model_parameters)
    output_replicate_path = joinpath(output_path, parameters.parameters, string(this_replicate_id))
    if ~ispath(output_replicate_path)
        mkpath(output_replicate_path)
    end

    # We ONLY store the raw output at this point, the processing of data is done in another job
    writedlm(joinpath(output_replicate_path, "timeseries.dat"), output[1])
    writedlm(joinpath(output_replicate_path, "mutualism.dat"), output[2])
    writedlm(joinpath(output_replicate_path, "competition.dat"), output[3])
    writedlm(joinpath(output_replicate_path, "predation.dat"), output[4])
    writedlm(joinpath(output_replicate_path, "infection.dat"), output[5])
end

# Read the folder and write the results for this simulation
result_path = joinpath(dirname(Base.active_project()), "results")
if ~ispath(result_path)
    mkpath(result_path)
end

results = [DataFrame() for thr in 1:Threads.nthreads()]

simulations = readdir(joinpath(output_path, parameters.parameters), join=true)
Threads.@threads for simulation in simulations
    try
        simid, repid = splitpath(simulation)[(end-1):end]
        timeseries = readdlm(joinpath(simulation, "timeseries.dat"))
        competition = readdlm(joinpath(simulation, "competition.dat"))
        infection = readdlm(joinpath(simulation, "infection.dat"))
        mutualism = readdlm(joinpath(simulation, "mutualism.dat"))
        predation = readdlm(joinpath(simulation, "predation.dat"))

        endpoint = timeseries[2:end, end]

        sp = round(Int, length(endpoint) / 2)

        Sₜ = endpoint[1:sp]
        Iₜ = endpoint[(sp+1):end]
        Hₜ = Sₜ + Iₜ

        keep = findall(Hₜ .> 0.1)

        Sₜ = Sₜ[keep]
        Iₜ = Iₜ[keep]
        Hₜ = Hₜ[keep]
        pₜ = Hₜ ./ sum(Hₜ)

        prevalence = Iₜ ./ Hₜ
        degree = vec(sum(infection[keep, keep] .> 0.0, dims=1))

        push!(results[Threads.threadid()], (
            parameters=simid,
            replicate=repid,
            richness=length(Hₜ),
            susceptibles=sum(Sₜ),
            infectious=sum(Iₜ),
            population=sum(Hₜ),
            correlation=sum(Iₜ) == 0 ? 0.0 : cor(degree, prevalence),
            prevalence=sum(Iₜ) / sum(Hₜ),
            diversity=sum(-pₜ .* log.(pₜ))
        ))
    catch err
        @info err
        continue
    end
end

CSV.write(joinpath(result_path, "$(parameters.parameters).csv"), vcat(results...))