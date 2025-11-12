#!/bin/bash

#SBATCH --partition=a3,a3mega
#SBATCH --output="run_calibration.txt"
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=1

julia --project=experiments/ClimaEarth -e 'using Pkg; Pkg.develop(;path="."); Pkg.instantiate(;verbose=true)'
julia --project=experiments/ClimaEarth -e 'using Pkg; \
    Pkg.add(Pkg.PackageSpec(;name="ClimaAnalysis", rev="main")); \
    Pkg.add(Pkg.PackageSpec(;name="ClimaCalibrate", rev="main")); \
'

julia --project=experiments/ClimaEarth experiments/calibration/run_calibration.jl

