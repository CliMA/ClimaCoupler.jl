import ClimaComms
ClimaComms.@import_required_backends
import ArgParse
import ClimaCoupler
import ClimaCoupler: Utilities
import ClimaAtmos as CA
import Random

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))
Random.seed!(1234)

function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        "--config_file"
        help = "A yaml file used to set the configuration of the coupled model"
        arg_type = String
        default = nothing
        "--job_id"
        help = "A unique string ID for this job"
        arg_type = String
        default = "climaatmos_test"
    end
    return s
end

# Read in atmos configuration file from command line
parse_commandline(s) = ArgParse.parse_args(ARGS, s)
parsed_args = parse_commandline(argparse_settings())
atmos_config_file = parsed_args["config_file"]
job_id = parsed_args["job_id"]

# Use input atmos configuration file
if isnothing(atmos_config_file)
    @info "Using Atmos default configuration"
    config = CA.default_config_dict()
else
    @info "Using Atmos configuration from $atmos_config_file"
    config = CA.override_default_config(joinpath(pkgdir(ClimaCoupler), atmos_config_file))
end

# Specify atmos output directory to be inside the coupler output directory
output_dir = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "ClimaEarth",
    "output",
    job_id,
    "artifacts",
)
!isdir(output_dir) && mkpath(output_dir)
config = merge(config, Dict("output_dir" => output_dir))
atmos_config = CA.AtmosConfig(config)
simulation = CA.get_simulation(atmos_config)
sol_res = CA.solve_atmos!(simulation)

## Use ClimaAtmos calculation to show the simulated years per day of the simulation (SYPD)
tspan = (Float64(0), Float64(Utilities.time_to_seconds(config["t_end"])))
walltime = sol_res.walltime
es = CA.EfficiencyStats(tspan, walltime)
sypd = CA.simulated_years_per_day(es)
@info "SYPD: $sypd"

## Save the SYPD and max RSS information
comms_ctx = atmos_config.comms_ctx
if ClimaComms.iamroot(comms_ctx)
    open(joinpath(output_dir, "sypd.txt"), "w") do sypd_filename
        println(sypd_filename, "$sypd")
    end

    open(joinpath(output_dir, "max_rss_cpu.txt"), "w") do cpu_max_rss_filename
        cpu_max_rss_GB = Utilities.show_memory_usage()
        println(cpu_max_rss_filename, cpu_max_rss_GB)
    end
end
