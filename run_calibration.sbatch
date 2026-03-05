#!/bin/bash
#SBATCH --job-name=slurm_calibration
#SBATCH --output=calibration_%j.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G

# Assuming the script is run from the top level of ClimaCoupler.jl
CLIMAEARTH_PATH='experiments/calibration/subseasonal/'
DRIVER=${CLIMAEARTH_PATH}'run_calibration.jl'

# Set environment variables for CliMA
# Note: CUDA device is set by ClimaCalibrate's spawned member jobs
export CLIMACOMMS_CONTEXT="SINGLETON"

# Load required modules
module load climacommon

# Build and run the Julia code
julia --project=calibration -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.develop(Pkg.PackageSpec(;name="ClimaAtmos", path="/home/chew/git_repos/ClimaAtmos.jl"))'
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="ClimaParams", rev="rc/gw"))'
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.precompile()'
julia --project=experiments/ClimaEarth $DRIVER
