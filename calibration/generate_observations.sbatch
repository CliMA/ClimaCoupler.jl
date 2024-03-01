#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=8
#SBATCH --output="experiments/amip_coupled/truth_simulation/model_log.out"

# Configure the environment
source experiments/common_env.sh

echo "Generating truth observations."

srun  --open-mode=append julia --project=experiments -e '
using ClimaComms
ClimaComms.init(ClimaComms.context())
import ClimaCoupler as CCo
import YAML
using NetCDF
import JLD2

experiment_dir = joinpath("experiments", "amip_coupled")
output_dir = joinpath(experiment_dir, "truth_simulation")

config_dict = YAML.load_file(joinpath(experiment_dir, "coupler_config.yml"))
config_dict["t_end"] = "300days"
config_dict["output_dir"] = output_dir


coupled_simulation = 
CA.solve_coupler!(coupled_simulation)

ts = ncread(joinpath(output_dir, "ts_60d_average.nc"), "ta")
include(joinpath(experiment_dir, "observation_map.jl"))
(; observation, variance) = process_member_data(ta; output_variance = true)
JLD2.save_object(joinpath(experiment_dir, "obs_mean.jld2"), observation)
JLD2.save_object(joinpath(experiment_dir, "obs_noise_cov.jld2"), variance)
'
