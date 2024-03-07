## coupler defaults
# get component model dictionaries
include("../experiments/AMIP/cli_options.jl")
parsed_args = parse_commandline(argparse_settings())
config_dict = YAML.load_file("./experiments/amip_coupled/coupler_config.yml")
config_dict = YAML.load_file(joinpath(experiment_dir, "coupler_config.yml"));
config_dict["t_end"] = "150secs";
config_dict["output_dir"] = output_dir;
config_dict = merge(parsed_args, config_dict)
config_dict_atmos = get_atmos_config(config_dict)

# merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
# (if there are common keys, the last dictorionary in the `merge` arguments takes precedence)
config_dict = merge(config_dict_atmos, config_dict)


## read in some parsed command line arguments
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
date0 = date = DateTime(config_dict["start_date"], dateformat"yyyymmdd")
mono_surface = config_dict["mono_surface"]
hourly_checkpoint = config_dict["hourly_checkpoint"]
restart_dir = config_dict["restart_dir"]
restart_t = Int(config_dict["restart_t"])
evolving_ocean = config_dict["evolving_ocean"]
config_dict["print_config_dict"] = false

## I/O directory setup
COUPLER_OUTPUT_DIR = "/Users/akshaysridhar/Research/Codes/ClimaCoupler.jl/calibration/output/amip/"
mkpath(COUPLER_OUTPUT_DIR)

REGRID_DIR = joinpath(COUPLER_OUTPUT_DIR, "regrid_tmp/")
mkpath(REGRID_DIR)

COUPLER_ARTIFACTS_DIR = COUPLER_OUTPUT_DIR * "_artifacts"
isdir(COUPLER_ARTIFACTS_DIR) ? nothing : mkpath(COUPLER_ARTIFACTS_DIR)

config_dict["print_config_dict"] ? @info(config_dict) : nothing
config_dict_atmos["output_dir"] = COUPLER_OUTPUT_DIR
