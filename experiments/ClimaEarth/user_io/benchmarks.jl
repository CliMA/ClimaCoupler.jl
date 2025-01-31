#=
Our goal here is to output a table displaying some results from benchmark runs
in the coupler. We want to be able to compare between CPU and GPU runs, as well
as between coupled and atmos-only runs. The metrics we want to compare are
SYPD, memory usage, allocations, and the maximum, median, and mean differences
between the CPU and GPU states.

The table should look something like this (note that the last 3 columns will be
added in a future PR):
------------------------------------------------------------------------------------
|             | CPU Run      | GPU Run     | Max Diff  |  Median Diff |  Mean Diff |
------------------------------------------------------------------------------------
|             |   $job_id    |   $job_id   |           |              |              |
| Coupled run |    $SYPD     |    $SYPD    | $max_diff | $median_diff |  $mean_diff  |
|             | $cpu_max_rss | $cpu_max_rss|           |              |              |
------------------------------------------------------------------------------------
|             |   $job_id    |   $job_id   |           |              |              |
| Atmos-only  |    $SYPD     |    $SYPD    | $max_diff | $median_diff |  $mean_diff  |
|             | $cpu_max_rss | $cpu_max_rss|           |              |              |
------------------------------------------------------------------------------------

=#

import ArgParse
import PrettyTables
import ClimaCoupler

function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        "--cpu_job_id_coupled"
        help = "The name of the CPU coupled run we want to compare. User must specify CPU and/or GPU coupled run name."
        arg_type = String
        default = nothing
        "--gpu_job_id_coupled"
        help = "The name of the GPU coupled run we want to compare."
        arg_type = String
        default = nothing
        "--cpu_job_id_coupled_io"
        help = "The name of the CPU coupled run with IO we want to compare. User must specify CPU and/or GPU coupled run name."
        arg_type = String
        default = nothing
        "--gpu_job_id_coupled_io"
        help = "The name of the GPU coupled with IO run we want to compare."
        arg_type = String
        default = nothing
        "--cpu_job_id_atmos"
        help = "The name of the CPU atmos-only run without diagnostic EDMF we want to compare. User must specify CPU and/or GPU atmos-only non-EDMF run name."
        arg_type = String
        default = nothing
        "--gpu_job_id_atmos"
        help = "The name of the GPU atmos-only run we want to compare."
        arg_type = String
        default = nothing
        "--cpu_job_id_atmos_diagedmf"
        help = "The name of the CPU atmos-only run with diagnostic EDMF we want to compare. User must specify CPU and/or GPU atmos-only EDMF run name."
        arg_type = String
        default = nothing
        "--gpu_job_id_atmos_diagedmf"
        help = "The name of the GPU atmos-only run we want to compare."
        arg_type = String
        default = nothing
        "--coupler_output_dir"
        help = "Directory to save output files."
        arg_type = String
        default = "experiments/ClimaEarth/output"
        "--build_id"
        help = "The build ID of the pipeline running this script."
        arg_type = String
        default = nothing
    end
    return s
end

"""
    get_run_info(parsed_args, run_type)

Use the input `parsed_args` to get the job ID and artifacts directories for
both the CPU and GPU runs of the given `run_type`.

`run_type` must be one of "coupled", "coupled_io", "atmos", or "atmos_diagedmf".
"""
function get_run_info(parsed_args, run_type)
    # Read in CPU and GPU job ID info from command line
    if run_type == "coupled"
        cpu_job_id = parsed_args["cpu_job_id_coupled"]
        gpu_job_id = parsed_args["gpu_job_id_coupled"]
        mode_name = "amip"
    elseif run_type == "coupled_io"
        cpu_job_id = parsed_args["cpu_job_id_coupled_io"]
        gpu_job_id = parsed_args["gpu_job_id_coupled_io"]
        mode_name = "amip"
    elseif run_type == "atmos_diagedmf"
        cpu_job_id = parsed_args["cpu_job_id_atmos_diagedmf"]
        gpu_job_id = parsed_args["gpu_job_id_atmos_diagedmf"]
        mode_name = "climaatmos"
    elseif run_type == "atmos"
        cpu_job_id = parsed_args["cpu_job_id_atmos"]
        gpu_job_id = parsed_args["gpu_job_id_atmos"]
        mode_name = "climaatmos"
    else
        error("Invalid run type: $run_type")
    end

    # Verify that the user has provided the necessary job IDs
    # If only one job ID of the CPU/GPU run pair is provided, the other will be inferred
    if isnothing(cpu_job_id) && isnothing(gpu_job_id)
        error("Must pass CPU and/or GPU coupled run name to compare them.")
    elseif isnothing(gpu_job_id)
        gpu_job_id = "gpu_" * cpu_job_id
    elseif isnothing(cpu_job_id)
        cpu_job_id = gpu_job_id[5:end]
    end

    # Construct CPU and GPU artifacts directories
    cpu_artifacts_dir = joinpath(output_dir, cpu_job_id, "artifacts")
    gpu_artifacts_dir = joinpath(output_dir, gpu_job_id, "artifacts")

    return (cpu_job_id, gpu_job_id, cpu_artifacts_dir, gpu_artifacts_dir)
end

"""
    get_run_data(artifacts_dir)

Read in run data from artifacts directories, currently SYPD and max RSS on the CPU.
"""
function get_run_data(artifacts_dir)
    # Read in SYPD info
    sypd = open(joinpath(artifacts_dir, "sypd.txt"), "r") do sypd_file
        round(parse(Float64, read(sypd_file, String)), digits = 4)
    end

    # Read in max RSS info
    cpu_max_rss = open(joinpath(artifacts_dir, "max_rss_cpu.txt"), "r") do cpu_max_rss_file
        read(cpu_max_rss_file, String)
    end

    return (sypd, cpu_max_rss)
end

"""
    append_table_data(table_data, setup_id, cpu_job_id, gpu_job_id, cpu_artifacts_dir, gpu_artifacts_dir)

Append data for a given setup to the table data.
"""
function append_table_data(table_data, setup_id, cpu_job_id, gpu_job_id, cpu_artifacts_dir, gpu_artifacts_dir)
    # Get SYPD and allocation info for both input runs
    cpu_sypd, cpu_max_rss = get_run_data(cpu_artifacts_dir)
    gpu_sypd, gpu_cpu_max_rss = get_run_data(gpu_artifacts_dir)

    # Create rows containing data for these runs
    new_table_data = [
        ["" "job ID:" cpu_job_id gpu_job_id]
        [setup_id "SYPD:" cpu_sypd gpu_sypd]
        ["" "CPU max RSS:" cpu_max_rss gpu_cpu_max_rss]
    ]
    return vcat(table_data, new_table_data)
end


# Read in command line arguments
parsed_args = ArgParse.parse_args(ARGS, argparse_settings())
output_dir = parsed_args["coupler_output_dir"]

# Access buildkite pipeline ID (from `BUILDKITE_GITHUB_DEPLOYMENT_ID` variable)
build_id = parsed_args["build_id"]
if !isnothing(build_id)
    build_id_str = "Build ID: $build_id"
else
    build_id_str = "Build ID: N/A"
end

# Read in run info for each of the cases we want to compare
run_info_coupled = get_run_info(parsed_args, "coupled")
run_info_coupled_io = get_run_info(parsed_args, "coupled_io")
run_info_atmos_diagedmf = get_run_info(parsed_args, "atmos_diagedmf")
run_info_atmos = get_run_info(parsed_args, "atmos")

# Set up info for PrettyTables.jl
headers = [build_id_str, "Horiz. res.: 30 elems", "CPU Run [64 processes]", "GPU Run [2 A100s]"]
data = [
    ["" "Vert. res.: 63 levels" "" ""]
    ["" "dt: 120secs" "" ""]
]

# Append data to the table for each of the cases we want to compare
data = append_table_data(data, "Coupled with diag. EDMF + IO", run_info_coupled_io...)
data = append_table_data(data, "Coupled with diag. EDMF", run_info_coupled...)
data = append_table_data(data, "Atmos with diag. EDMF", run_info_atmos_diagedmf...)
data = append_table_data(data, "Atmos without diag. EDMF", run_info_atmos...)

# Use the coupled CPU job ID for the output dir
cpu_job_id_coupled = run_info_coupled[1]
table_output_dir = joinpath(output_dir, "compare_amip_climaatmos_$(cpu_job_id_coupled)")
!isdir(table_output_dir) && mkdir(table_output_dir)
table_path = joinpath(table_output_dir, "table.txt")
open(table_path, "w") do f
    # Output the table, including lines before and after the header
    PrettyTables.pretty_table(f, data, header = headers, hlines = [0, 3, 6, 9, 12, 15])
end
