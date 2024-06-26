### Generate synthetic truth datasets
@info "Generating synthetic truth data"
using ClimaComms
import ClimaCalibrate
ClimaComms.init(ClimaComms.context())
import ClimaCoupler 
import YAML
using NCDatasets
import JLD2
using Statistics

import ClimaCalibrate: get_forward_model

experiment_dir = joinpath("experiments", "amip_coupled")
COUPLER_OUTPUT_DIR = joinpath(experiment_dir, "truth_simulation")
include("coupler_driver_calibration.jl");
cs = get_simulation(config_dict);
solve_coupler!(cs); # Integrate the coupled model

### Process "Observations" -> Store in `testdir`
testdir = "/Users/akshaysridhar/Research/Codes/ClimaCoupler.jl/calibration/experiments/amip_coupled/truth_simulation/"
wa = NCDataset(joinpath(testdir, "", "wa_inst.nc"))["wa"]
include(joinpath(experiment_dir, "observation_map.jl"))
(; observation, variance) = process_member_data(wa; output_variance = true)
JLD2.save_object(joinpath(experiment_dir, "obs_mean.jld2"), observation)
JLD2.save_object(joinpath(experiment_dir, "obs_noise_cov.jld2"), variance)

### Run forward model iterations
include("coupler_interface.jl")
experiment_id = "amip_coupled"

iteration = 1
format_i = "iteration_$iteration"

SLURM_ARRAY_TASK_ID = 1
member = "member_$SLURM_ARRAY_TASK_ID"
output="output/$experiment_id/$format_i/$member/model_log.out"

experiment_dir = joinpath("experiments", "amip_coupled")
COUPLER_OUTPUT_DIR = joinpath("experiments","$format_i","$member")
(;config_dict) = get_config(CoupledModel(), 1, iteration, experiment_id);
run_forward_model(CoupledModel(), config_dict); 
# This outputs in experiments/AMIP/output/amip/target_amip_n1_shortrun

### Calibrate
ClimaCalibrate.calibrate("amip_coupled")

### Re-run target simulation
