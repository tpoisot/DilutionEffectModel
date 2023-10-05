#! /bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2500M
#SBATCH --cpus-per-task=64
#SBATCH --job-name=disease-dilution-dataprep
#SBATCH --output=slurm/%x-%a.out

module load StdEnv/2020 julia/1.9.1
julia --project -t 64 02_read_simulations.jl
