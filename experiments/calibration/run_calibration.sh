#!/bin/bash

#SBATCH --partition=a3
#SBATCH --output="run_calibration.txt"
#SBATCH --time=05:00:00
#SBATCH --ntasks=10
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=4

julia --project=experiments/AMIP -e 'using Pkg; Pkg.develop(;path="."); Pkg.instantiate(;verbose=true)'

julia --project=experiments/AMIP experiments/calibration/run_calibration.jl
