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
parameters = parameter_sets[parameter_set_id, :]

# Seed the RNG for the node
Random.seed!(parameters.seed)

# Prepare temporary data frames for each thread
results = [DataFrame() for thr in 1:Threads.nthreads()]

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

Threads.@threads for replicate in 1:parameters.replicates
    this_replicate_id = uuid4()
    output = onesim(parameters.richness, interaction_parameters, model_parameters)
    output_replicate_path = joinpath(output_path, parameters.parameters, string(this_replicate_id))
    if ~ispath(output_replicate_path)
        mkpath(output_replicate_path)
    end
    writedlm(joinpath(output_replicate_path, "timeseries.dat"), output[1])
    writedlm(joinpath(output_replicate_path, "mutualism.dat"), output[2])
    writedlm(joinpath(output_replicate_path, "competition.dat"), output[3])
    writedlm(joinpath(output_replicate_path, "predation.dat"), output[4])
    writedlm(joinpath(output_replicate_path, "infection.dat"), output[5])
end