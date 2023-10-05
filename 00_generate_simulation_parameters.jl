using UUIDs
using DataFrames
import CSV

input_path = joinpath(dirname(Base.active_project()), "inputs")
if ~ispath(input_path)
    mkpath(input_path)
end

output_path = joinpath(dirname(Base.active_project()), "outputs")
if ~ispath(output_path)
    mkpath(output_path)
end

if ~ispath("slurm")
    mkpath("slurm")
end

bank = DataFrame()
centers = collect(0:2:100)
barycenters = [(c1, c2, c3) for c1 in centers for c2 in centers for c3 in centers]
filter!(b -> isequal(100)(sum(b)), barycenters)

for b in barycenters
    for density in [true, false]
        sigma_infection = density ? 0.001 : 0.01
        push!(bank, (
            parameters=uuid4(),
            seed=rand(1:200_000),
            mutualism=b[1],
            competition=b[2],
            predation=b[3],
            density=density,
            replicates=1000,
            richness=50,
            time=1000,
            connectance=0.5,
            carrying_capacity=100.0,
            sigma_interaction=0.2,
            sigma_infection=sigma_infection,
            regulation=0.02,
            virulence=0.05,
            recovery=0.05,
            half_saturation=30.0
        ))
    end
end

CSV.write(joinpath(input_path, "sets.csv"), bank)

# Write the SLURM file to go with the job
job_file = """
#! /bin/bash
#SBATCH --array=1-$(size(bank, 1))
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2500M
#SBATCH --cpus-per-task=64
#SBATCH --job-name=disease-dilution-effect
#SBATCH --output=$(joinpath("slurm", "%x-%a.out"))

module load StdEnv/2020 julia/1.9.1
julia --project -t 64 01_run_one_parameter.jl
"""
write("dilution.sh", job_file)