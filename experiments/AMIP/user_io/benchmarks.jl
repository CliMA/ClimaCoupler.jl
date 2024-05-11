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
|             | $cpu_allocs  | $cpu_allocs |           |              |              |
------------------------------------------------------------------------------------
|             |  $run_name   |  $run_name  |           |              |              |
| Atmos-only  |    $SYPD     |    $SYPD    | $max_diff | $median_diff |  $mean_diff  |
|             | $cpu_allocs  | $cpu_allocs |           |              |              |
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
        "--build_id"
        help = "The build ID of the pipeline running this script."
        arg_type = String
        default = nothing
    end
    return s
end

# Parse command line arguments
parsed_args = ArgParse.parse_args(ARGS, argparse_settings())

# Access buildkite pipeline ID (from `BUILDKITE_GITHUB_DEPLOYMENT_ID` variable)
build_id = parsed_args["build_id"]
if !isnothing(build_id)
    build_id_str = "Build ID: $build_id"
else
    build_id_str = "Build ID: N/A"
end

# Construct CPU and GPU artifacts directories
output_dir = parsed_args["coupler_output_dir"]

# Coupled runs
# Read in CPU and GPU run name info from command line
cpu_run_name_coupled = parsed_args["cpu_run_name_coupled"]
gpu_run_name_coupled = parsed_args["gpu_run_name_coupled"]
if isnothing(cpu_run_name_coupled) && isnothing(gpu_run_name_coupled)
    error("Must pass CPU and/or GPU coupled run name to compare them.")
elseif isnothing(gpu_run_name_coupled)
    gpu_run_name_coupled = "gpu_" * cpu_run_name_coupled
elseif isnothing(cpu_run_name_coupled)
    cpu_run_name_coupled = gpu_run_name_coupled[5:end]
end

# Read in mode name from command line (or retrieve from run name).
# Note that we expect this to be the same for all 4 simulations being compared.
mode_name = parsed_args["mode_name"]
if isnothing(mode_name)
    mode_name =
        occursin("amip", cpu_run_name_coupled) ? "amip" :
        (occursin("slabplanet", cpu_run_name_coupled) ? "slabplanet" : error("Please provide a valid `mode_name`."))
end

gpu_artifacts_dir_coupled = joinpath(output_dir, mode_name, gpu_run_name_coupled) * "_artifacts"
cpu_artifacts_dir_coupled = joinpath(output_dir, mode_name, cpu_run_name_coupled) * "_artifacts"

# Atmos-only runs
# Read in CPU and GPU run name info from command line
cpu_run_name_atmos = parsed_args["cpu_run_name_atmos"]
gpu_run_name_atmos = parsed_args["gpu_run_name_atmos"]
if isnothing(cpu_run_name_atmos) && isnothing(gpu_run_name_atmos)
    error("Must pass CPU and/or GPU coupled run name to compare them.")
elseif isnothing(gpu_run_name_atmos)
    gpu_run_name_atmos = "gpu_" * cpu_run_name_atmos
elseif isnothing(cpu_run_name_atmos)
    cpu_run_name_atmos = gpu_run_name_atmos[5:end]
    cpu_artifacts_dir_atmos = joinpath(output_dir, cpu_run_name_atmos)
end

mode_name_atmos = "climaatmos"
gpu_artifacts_dir_atmos = joinpath(output_dir, mode_name_atmos, gpu_run_name_atmos) * "_artifacts"
cpu_artifacts_dir_atmos = joinpath(output_dir, mode_name_atmos, cpu_run_name_atmos) * "_artifacts"

# Read in SYPD and allocations info from artifacts directories
function get_sypd_allocs(artifacts_dir)
    # Read in SYPD info
    sypd_file = open(joinpath(artifacts_dir, "sypd.txt"), "r")
    sypd = round(parse(Float64, read(sypd_file, String)), digits = 4)

    # Read in allocations info
    cpu_allocs_file = open(joinpath(artifacts_dir, "allocations_cpu.txt"), "r")
    cpu_allocs = read(cpu_allocs_file, String)

    return (sypd, cpu_allocs)
end

cpu_sypd_coupled, cpu_allocs_coupled = get_sypd_allocs(cpu_artifacts_dir_coupled)
gpu_sypd_coupled, gpu_cpu_allocs_coupled = get_sypd_allocs(gpu_artifacts_dir_coupled)
cpu_sypd_atmos, cpu_allocs_atmos = get_sypd_allocs(cpu_artifacts_dir_atmos)
gpu_sypd_atmos, gpu_cpu_allocs_atmos = get_sypd_allocs(gpu_artifacts_dir_atmos)

# Set up info for PrettyTables.jl
headers = [build_id_str, "Horiz. res.: 30 elems", "CPU Run [64 processes]", "GPU Run [4 A100s]"]
data = [
    ["" "Vert. res.: 63 levels" "" ""]
    ["" "dt: 120secs" "" ""]
    ["" "run name:" cpu_run_name_coupled gpu_run_name_coupled]
    ["Coupled" "SYPD:" cpu_sypd_coupled gpu_sypd_coupled]
    ["" "CPU max RSS allocs:" cpu_allocs_coupled gpu_cpu_allocs_coupled]
    ["" "run name:" cpu_run_name_atmos gpu_run_name_atmos]
    ["Atmos-only" "SYPD:" cpu_sypd_atmos gpu_sypd_atmos]
    ["" "CPU max RSS allocs:" cpu_allocs_atmos gpu_cpu_allocs_atmos]
]

# Use the coupled CPU run name for the output dir
table_output_dir = joinpath(output_dir, "compare_$(mode_name)_$(mode_name_atmos)_$(cpu_run_name_coupled)")
!isdir(table_output_dir) && mkdir(table_output_dir)
table_path = joinpath(table_output_dir, "table.txt")
open(table_path, "w") do f
    # Output the table, including lines before and after the header
    PrettyTables.pretty_table(f, data, header = headers, hlines = [0, 3, 6, 9])
end
