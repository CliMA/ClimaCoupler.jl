#!/bin/bash

#SBATCH --partition=a3
#SBATCH --output="run_calibration.txt"
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=1

julia --project=experiments/ClimaEarth -e 'using Pkg; Pkg.develop(;path="."); Pkg.instantiate(;verbose=true)'

julia --project=experiments/ClimaEarth experiments/calibration/coarse_amip/run_calibration.jl

