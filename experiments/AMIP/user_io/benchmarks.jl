using PrettyTables

#=
Our goal here is to output a table displaying some results from benchmark runs
in the coupler. We want to be able to compare between CPU and GPU runs, as well
as between coupled and atmos-only runs. The metrics we want to compare are
SYPD, allocations, and the maximum, median, and mean differences between the
CPU and GPU states.

The table should look something like this (note that the last 3 columns will be
added in a future PR):
------------------------------------------------------------------------------------
|             | CPU Run      | GPU Run     | Max Diff  |  Median Diff |  Mean Diff |
------------------------------------------------------------------------------------
|             |  $run_name   |  $run_name  |           |              |              |
| Coupled run |    $SYPD     |    $SYPD    | $max_diff | $median_diff |  $mean_diff  |
|             | $cpu_allocs  | $gpu_allocs |           |              |              |
------------------------------------------------------------------------------------
|             |  $run_name   |  $run_name  |           |              |              |
| Atmos-only  |    $SYPD     |    $SYPD    | $max_diff | $median_diff |  $mean_diff  |
|             | $cpu_allocs  | $gpu_allocs |           |              |              |
------------------------------------------------------------------------------------

=#

import ArgParse
import PrettyTables
import ClimaCoupler

function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        "--cpu_run_name_coupled"
        help = "The name of the CPU coupled run we want to compare. User must specify CPU and/or GPU coupled run name."
        arg_type = String
        default = nothing
        "--gpu_run_name_coupled"
        help = "The name of the GPU coupled run we want to compare."
        arg_type = String
        default = nothing
        "--cpu_run_name_atmos"
        help = "The name of the CPU atmos-only run we want to compare. User must specify CPU and/or GPU atmos-only run name."
        arg_type = String
        default = nothing
        "--gpu_run_name_atmos"
        help = "The name of the GPU atmos-only run we want to compare."
        arg_type = String
        default = nothing
        "--mode_name"
        help = "The mode of the simulations being compared (`slabplanet` or `AMIP`)."
        arg_type = String
        default = nothing
        "--coupler_output_dir"
        help = "Directory to save output files."
        arg_type = String
        default = "experiments/AMIP/output"
    end
    return s
end

# Read in CPU and GPU run name info from command line
parsed_args = ArgParse.parse_args(ARGS, argparse_settings())

# Coupled runs
cpu_run_name_coupled = parsed_args["cpu_run_name_coupled"]
gpu_run_name_coupled = parsed_args["gpu_run_name_coupled"]
if isnothing(cpu_run_name_coupled) && isnothing(gpu_run_name_coupled)
    error("Must pass CPU and/or GPU coupled run name to compare them.")
elseif isnothing(gpu_run_name_coupled)
    gpu_run_name_coupled = "gpu_" * cpu_run_name_coupled
elseif isnothing(cpu_run_name_coupled)
    cpu_run_name_coupled = gpu_run_name_coupled[5:end]
end

# Atmos-only runs
cpu_run_name_atmos = parsed_args["cpu_run_name_atmos"]
gpu_run_name_atmos = parsed_args["gpu_run_name_atmos"]
if isnothing(cpu_run_name_atmos) && isnothing(gpu_run_name_atmos)
    error("Must pass CPU and/or GPU coupled run name to compare them.")
elseif isnothing(gpu_run_name_atmos)
    gpu_run_name_atmos = "gpu_" * cpu_run_name_atmos
elseif isnothing(cpu_run_name_atmos)
    cpu_run_name_atmos = gpu_run_name_atmos[5:end]
end

# Read in mode name from command line (or retrieve from run name).
#Note that we expect this to be the same for all 4 simulations being compared.
mode_name = parsed_args["mode_name"]
if isnothing(mode_name)
    mode_name =
        occursin("amip", cpu_run_name) ? "amip" :
        (occursin("slabplanet", cpu_run_name) ? "slabplanet" : error("Please provide a valid `mode_name`."))
end

# Construct CPU and GPU artifacts directories
output_dir = parsed_args["coupler_output_dir"]

# Read in SYPD and allocations info from artifacts directories
function get_sypd_allocs(output_dir, cpu_run_name, gpu_run_name)
    cpu_artifacts_dir = joinpath(output_dir, joinpath(mode_name, cpu_run_name)) * "_artifacts"
    gpu_artifacts_dir = joinpath(output_dir, joinpath(mode_name, gpu_run_name)) * "_artifacts"

    # Read in SYPD info
    cpu_sypd_file = open(joinpath(cpu_artifacts_dir, "sypd.txt"), "r")
    cpu_sypd = read(cpu_sypd_file, String)

    gpu_sypd_file = open(joinpath(gpu_artifacts_dir, "sypd.txt"), "r")
    gpu_sypd = read(gpu_sypd_file, String)

    # Read in allocations info
    cpu_allocs_file = open(joinpath(cpu_artifacts_dir, "allocations.txt"), "r")
    cpu_allocs = read(cpu_allocs_file, String)

    gpu_allocs_file = open(joinpath(gpu_artifacts_dir, "allocations.txt"), "r")
    gpu_allocs = read(gpu_allocs_file, String)

    return (cpu_sypd, gpu_sypd, cpu_allocs, gpu_allocs)
end

cpu_sypd_coupled, gpu_sypd_coupled, cpu_allocs_coupled, gpu_allocs_coupled =
    get_sypd_allocs(output_dir, cpu_run_name_coupled, gpu_run_name_coupled)
cpu_sypd_atmos, gpu_sypd_atmos, cpu_allocs_atmos, gpu_allocs_atmos =
    get_sypd_allocs(output_dir, cpu_run_name_atmos, gpu_run_name_atmos)

# Set up info for PrettyTables.jl
headers = ["", "CPU Run", "GPU Run"]
data = [
    ["", cpu_run_name_coupled, gpu_run_name_coupled],
    ["Coupled run", cpu_sypd_coupled, gpu_sypd_coupled],
    ["", cpu_allocs_coupled, gpu_allocs_coupled],
    ["", cpu_run_name_atmos, gpu_run_name_atmos],
    ["Atmos-only", cpu_sypd_atmos, gpu_sypd_atmos],
    ["", cpu_allocs_atmos, gpu_allocs_atmos],
]

# Output the table, including lines before and after the header, and after the 4th row
PrettyTables.pretty_table(data, headers, hlines = [0, 1, 4])
