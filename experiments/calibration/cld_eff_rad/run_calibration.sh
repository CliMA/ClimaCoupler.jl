#!/bin/bash

#SBATCH --partition=a3
#SBATCH --output="run_calibration.txt"
#SBATCH --time=24:00:00
#SBATCH --ntasks=30
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=4

julia --project=experiments/calibration -e 'using Pkg; Pkg.develop(;path="."); Pkg.instantiate(;verbose=true)'

julia --project=experiments/calibration experiments/calibration/cld_eff_rad/run_calibration.jl

