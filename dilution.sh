#! /bin/bash
#SBATCH --array=1-2652
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2500M
#SBATCH --cpus-per-task=64
#SBATCH --job-name=disease-dilution-effect
#SBATCH --output=slurm/%x-%a.out

module load StdEnv/2020 julia/1.9.1
julia --project -t 64 01_run_one_parameter.jl
