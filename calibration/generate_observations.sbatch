#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=8
#SBATCH --output="experiments/amip_coupled/truth_simulation/model_log.out"

# Configure the environment
export MODULEPATH=/groups/esm/modules:$MODULEPATH
module load climacommon/2024_02_27

echo "Generating truth observations."

srun  --open-mode=append julia --project=experiments -e '
using ClimaComms
ClimaComms.init(ClimaComms.context())
import ClimaCoupler as CCo
import YAML
using NCDatasets
import JLD2
using Statistics

experiment_dir = joinpath("experiments", "amip_coupled")
COUPLER_OUTPUT_DIR = joinpath(experiment_dir, "truth_simulation")

include("coupler_driver_calibration.jl");
solve_coupler!(cs); # Integrate the coupled model

testdir = "/Users/akshaysridhar/Research/Codes/ClimaCoupler.jl/calibration/experiments/amip_coupled/truth_simulation/"
wa = NCDataset(joinpath(testdir, "", "wa_inst.nc"))["wa"]
include(joinpath(experiment_dir, "observation_map.jl"))
(; observation, variance) = process_member_data(wa; output_variance = true)
JLD2.save_object(joinpath(experiment_dir, "obs_mean.jld2"), observation)
JLD2.save_object(joinpath(experiment_dir, "obs_noise_cov.jld2"), variance)
'
