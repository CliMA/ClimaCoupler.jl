#=
Our goal here is to output a table displaying some results from benchmark runs
in the coupler. We mainly want to compare performance (via SYPD) between
coupled and atmos-only runs.

The table should look something like this:
-------------------------------------------------------
| Build ID: XYZ | Horiz. res: ___ | GPU Run [X A100s] |
|               | Vert. res: ___  |                   |
|               | dt: ______      |                   |
-------------------------------------------------------
|               |   job ID:       |   $job_id         |
| Coupled run   |   SYPD:         |    $SYPD          |
-------------------------------------------------------
|               |   job ID:       |   $job_id         |
| Atmos-only    |   SYPD:         |    $SYPD          |
-------------------------------------------------------
=#

import ArgParse
import PrettyTables
import ClimaCoupler

function benchmarks_argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        "--gpu_job_id_coupled"
        help = "The name of the GPU coupled run we want to compare."
        arg_type = String
        default = nothing
        "--gpu_job_id_coupled_io"
        help = "The name of the GPU coupled with IO run we want to compare."
        arg_type = String
        default = nothing
        "--gpu_job_id_atmos"
        help = "The name of the GPU atmos-only run we want to compare."
        arg_type = String
        default = nothing
        "--gpu_job_id_atmos_diagedmf"
        help = "The name of the GPU atmos-only run we want to compare."
        arg_type = String
        default = nothing
        "--gpu_job_id_coupled_progedmf_coarse"
        help = "The name of the coarser GPU coupled with prognostic EDMF + 1M run we want to compare."
        arg_type = String
        default = nothing
        "--gpu_job_id_coupled_progedmf_fine"
        help = "The name of the finer GPU coupled with prognostic EDMF + 1M run we want to compare."
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
the GPU run of the given `run_type`.

`run_type` must be one of "coupled", "coupled_io", "atmos", or "atmos_diagedmf".
"""
function get_run_info(parsed_args, run_type)
    # Read in GPU job ID info from command line
    if run_type == "coupled"
        gpu_job_id = parsed_args["gpu_job_id_coupled"]
        mode_name = "amip"
    elseif run_type == "coupled_io"
        gpu_job_id = parsed_args["gpu_job_id_coupled_io"]
        mode_name = "amip"
    elseif run_type == "atmos_diagedmf"
        gpu_job_id = parsed_args["gpu_job_id_atmos_diagedmf"]
        mode_name = "climaatmos"
    elseif run_type == "atmos"
        gpu_job_id = parsed_args["gpu_job_id_atmos"]
        mode_name = "climaatmos"
    elseif run_type == "coupled_progedmf_coarse"
        gpu_job_id = parsed_args["gpu_job_id_coupled_progedmf_coarse"]
        mode_name = "amip"
    elseif run_type == "coupled_progedmf_fine"
        gpu_job_id = parsed_args["gpu_job_id_coupled_progedmf_fine"]
        mode_name = "amip"
    else
        error("Invalid run type: $run_type")
    end

    # Verify that the user has provided the necessary job IDs
    isnothing(gpu_job_id) && error("Must pass GPU coupled run name to compare it.")

    # Construct GPU artifacts directory
    gpu_artifacts_dir = joinpath(output_dir, gpu_job_id, "artifacts")

    return (gpu_job_id, gpu_artifacts_dir)
end

"""
    get_run_data(artifacts_dir)

Read in run data from artifacts directory, currently SYPD.
"""
function get_run_data(artifacts_dir)
    # Read in SYPD info
    sypd = open(joinpath(artifacts_dir, "sypd.txt"), "r") do sypd_file
        round(parse(Float64, read(sypd_file, String)), digits = 4)
    end

    return sypd
end

"""
    append_table_data(table_data, setup_id, gpu_job_id, gpu_artifacts_dir)

Append data for a given setup to the table data.
"""
function append_table_data(table_data, setup_id, gpu_job_id, gpu_artifacts_dir)
    # Get SYPD info for the GPU run
    gpu_sypd = get_run_data(gpu_artifacts_dir)

    # Create rows containing data for these runs
    new_table_data = [
        [setup_id "job ID:" gpu_job_id]
        ["" "SYPD:" gpu_sypd]
    ]
    return vcat(table_data, new_table_data)
end
# Read in command line arguments
parsed_args = ArgParse.parse_args(ARGS, benchmarks_argparse_settings())
output_dir = parsed_args["coupler_output_dir"]

# Access buildkite pipeline ID (from `BUILDKITE_GITHUB_DEPLOYMENT_ID` variable)
build_id = parsed_args["build_id"]
if !isnothing(build_id)
    build_id_str = "Build ID: $build_id"
else
    build_id_str = "Build ID: N/A"
end

# Read in run info for each of the cases we want to compare
run_info_coupled_progedmf_coarse = get_run_info(parsed_args, "coupled_progedmf_coarse")
run_info_coupled_progedmf_fine = get_run_info(parsed_args, "coupled_progedmf_fine")
run_info_coupled_io = get_run_info(parsed_args, "coupled_io")
run_info_coupled = get_run_info(parsed_args, "coupled")
run_info_atmos_diagedmf = get_run_info(parsed_args, "atmos_diagedmf")
run_info_atmos = get_run_info(parsed_args, "atmos")

# Set up info for PrettyTables.jl
headers = [build_id_str, "Horiz. res.: 30 elems", "GPU Run [2 A100s]"]
data = [
    ["" "Vert. res.: 63 levels" ""]
    ["" "dt: 120secs" ""]
]

# Append data to the table for each of the cases we want to compare
data = append_table_data(
    data,
    "Coupled with progedmf + 1M (16 helem)",
    run_info_coupled_progedmf_coarse...,
)
data = append_table_data(
    data,
    "Coupled with progedmf + 1M (30 helem)",
    run_info_coupled_progedmf_fine...,
)
data = append_table_data(data, "Coupled with diag. EDMF + IO", run_info_coupled_io...)
data = append_table_data(data, "Coupled with diag. EDMF", run_info_coupled...)
data = append_table_data(data, "Atmos with diag. EDMF", run_info_atmos_diagedmf...)
data = append_table_data(data, "Atmos without diag. EDMF", run_info_atmos...)

# Output table to file, note that this must match the slack upload command in the pipeline.yml file
table_output_dir = joinpath(output_dir, "compare_amip_climaatmos_amip_diagedmf")
mkpath(table_output_dir)
table_path = joinpath(table_output_dir, "table.txt")
open(table_path, "w") do f
    # Output the table, including lines before and after the header
    PrettyTables.pretty_table(
        f,
        data,
        header = headers,
        hlines = [0, 3, 5, 7, 9, 11, 13, 15],
    )
end
