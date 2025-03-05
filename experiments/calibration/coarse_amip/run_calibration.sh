#!/bin/bash

#SBATCH --partition=a3
#SBATCH --output="run_calibration.txt"
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=1
# Trap EXIT , should we trap SIGTERM as well?
# trap 'handle_exit' EXIT

# handle_exit() {
# if [ ! -f "output/surface_fluxes_perfect_model/iteration_004/eki_file.jld2" ]; then
# 	echo "Resubmitting due to incomplete calibration..."
# 	sbatch "$0"
# else
# 	echo "Calibration complete. No resubmission needed."
# fi
# }
julia --project=experiments/ClimaEarth -e 'using Pkg; Pkg.develop(;path="."); Pkg.instantiate(;verbose=true)'

julia --project=experiments/ClimaEarth experiments/calibration/coarse_amip/run_calibration.jl

