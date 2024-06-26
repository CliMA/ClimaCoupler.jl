#=
### Configuration Dictionaries
Each simulation mode has its own configuration dictionary. The `config_dict` of each simulation is a merge of the default configuration
dictionary and the simulation-specific configuration dictionary, which allows the user to override the default settings.

We can additionally pass the configuration dictionary to the component model initializers, which will then override the default settings of the component models.
=#

## coupler simulation default configuration
include("../experiments/AMIP/cli_options.jl")
parsed_args = parse_commandline(argparse_settings())

## read in config dictionary from file, overriding the coupler defaults
config_dict = YAML.load_file(joinpath(experiment_dir, "model_config.yml"));
config_dict = merge(parsed_args, config_dict)

## get component model dictionaries (if applicable)
config_dict_atmos = get_atmos_config(config_dict)

## merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
## (if there are common keys, the last dictorionary in the `merge` arguments takes precedence)
config_dict = merge(config_dict_atmos, config_dict)

## read in some parsed command line arguments, required by this script
mode_name = config_dict["mode_name"]
run_name = config_dict["run_name"]
energy_check = config_dict["energy_check"]
FT = config_dict["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim_name = "bucket"
t_end = Float64(time_to_seconds(config_dict["t_end"]))
t_start = 0.0
tspan = (t_start, t_end)
Δt_cpl = Float64(config_dict["dt_cpl"])
saveat = Float64(time_to_seconds(config_dict["dt_save_to_sol"]))
date0 = date = DateTime(config_dict["start_date"], "yyyymmdd")
mono_surface = config_dict["mono_surface"]
hourly_checkpoint = config_dict["hourly_checkpoint"]
restart_dir = config_dict["restart_dir"]
restart_t = Int(config_dict["restart_t"])
evolving_ocean = config_dict["evolving_ocean"]

#=
## Setup Communication Context
We set up communication context for CPU single thread/CPU with MPI/GPU. If no device is passed to `ClimaComms.context()`
then `ClimaComms` automatically selects the device from which this code is called.
=#

using ClimaComms
comms_ctx = ClimaCoupler.Utilities.get_comms_context(parsed_args)
const pid, nprocs = ClimaComms.init(comms_ctx)

#=
### I/O Directory Setup
`COUPLER_OUTPUT_DIR` is the directory where the output of the simulation will be saved, and `COUPLER_ARTIFACTS_DIR` is the directory where
the plots (from postprocessing and the conservation checks) of the simulation will be saved. `REGRID_DIR` is the directory where the regridding
temporary files will be saved.
=#

mkpath(COUPLER_OUTPUT_DIR)

REGRID_DIR = joinpath(COUPLER_OUTPUT_DIR, "regrid_tmp/")
mkpath(REGRID_DIR)

COUPLER_ARTIFACTS_DIR = COUPLER_OUTPUT_DIR * "_artifacts"
isdir(COUPLER_ARTIFACTS_DIR) ? nothing : mkpath(COUPLER_ARTIFACTS_DIR)

dir_paths = (; output = COUPLER_OUTPUT_DIR, artifacts = COUPLER_ARTIFACTS_DIR)

if ClimaComms.iamroot(comms_ctx)
    @info(COUPLER_OUTPUT_DIR)
    config_dict["print_config_dict"] ? @info(config_dict) : nothing
end