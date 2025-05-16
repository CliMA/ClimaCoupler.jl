#!/bin/bash

#SBATCH --partition=a3
#SBATCH --output="run_calibration.txt"
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1

julia --project=experiments/calibration -e 'using Pkg; Pkg.develop(;path="."); Pkg.instantiate(;verbose=true)'

julia --project=experiments/calibration experiments/calibration/coarse_amip/run_calibration.jl

