using ClimaComms
using Logging
import ArgParse

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))
import ClimaCoupler
import ClimaAtmos as CA
import Random
Random.seed!(1234)

function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        "--config_file"
        help = "A yaml file used to set the configuration of the coupled model"
        arg_type = String
        default = nothing
    end
    return s
end

# Read in atmos configuration file from command line
parsed_args = parse_commandline(argparse_settings())
atmos_config_file = parsed_args["config_file"]

# Use input atmos configuration file
if !(isnothing(atmos_config_file))
    @info "Using Atmos default configuration"
    config = CA.default_config_dict()
else
    @info "Using Atmos configuration from $atmos_config_file"
    config = CA.override_default_config(joinpath(pkgdir(CA), atmos_config_file))
end

# Specify atmos output directory to be inside the coupler output directory
output_dir = joinpath(pkgdir(ClimaCoupler), "experiments/AMIP/output", config["job_id"])
config = merge(config, Dict("output_dir" => output_dir))

simulation = CA.get_simulation(config)
(; integrator) = simulation
sol_res = CA.solve_atmos!(simulation)
