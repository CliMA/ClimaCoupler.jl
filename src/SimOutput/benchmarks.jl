#=
Our goal here is to output a table displaying some results from benchmark runs
in the coupler. We mainly want to compare SYPD between coupled and atmos-only runs.

The table should look something like this (note that the last 3 columns will be
added in a future PR):
-----------------------------
|             | GPU Run     |
-----------------------------
|             |   $job_id   |
| Coupled run |    $SYPD    |
-----------------------------
|             |   $job_id   |
| Atmos-only  |    $SYPD    |
-----------------------------

=#

import ArgParse
import ClimaCoupler

export get_benchmark_args, get_run_info, append_table_data

function argparse_settings_benchmarks()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        "--job_id_coupled"
        help = "The name of the GPU coupled run we want to compare."
        arg_type = String
        default = nothing
        "--job_id_coupled_io"
        help = "The name of the GPU coupled with IO run we want to compare."
        arg_type = String
        default = nothing
        "--job_id_atmos"
        help = "The name of the GPU atmos-only run we want to compare."
        arg_type = String
        default = nothing
        "--job_id_atmos_diagedmf"
        help = "The name of the GPU atmos-only run we want to compare."
        arg_type = String
        default = nothing
        "--job_id_coupled_progedmf_coarse"
        help = "The name of the coarser GPU coupled with prognostic EDMF + 1M run we want to compare."
        arg_type = String
        default = nothing
        "--job_id_coupled_progedmf_fine"
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
get_benchmark_args() = ArgParse.parse_args(ARGS, argparse_settings_benchmarks())

"""
    get_run_info(parsed_args, run_type)

Use the input `parsed_args` to get the job ID and artifacts directories for
the GPU run of the given `run_type`.

`run_type` must be one of "coupled", "coupled_io", "atmos", "atmos_diagedmf",
"coupled_progedmf_coarse", or "coupled_progedmf_fine".
"""
function get_run_info(parsed_args, run_type)
    # Read in GPU job ID info from command line
    if run_type == "coupled"
        job_id = parsed_args["job_id_coupled"]
        mode_name = "amip"
    elseif run_type == "coupled_io"
        job_id = parsed_args["job_id_coupled_io"]
        mode_name = "amip"
    elseif run_type == "atmos_diagedmf"
        job_id = parsed_args["job_id_atmos_diagedmf"]
        mode_name = "climaatmos"
    elseif run_type == "atmos"
        job_id = parsed_args["job_id_atmos"]
        mode_name = "climaatmos"
    elseif run_type == "coupled_progedmf_coarse"
        job_id = parsed_args["job_id_coupled_progedmf_coarse"]
        mode_name = "amip"
    elseif run_type == "coupled_progedmf_fine"
        job_id = parsed_args["job_id_coupled_progedmf_fine"]
        mode_name = "amip"
    else
        error("Invalid run type: $run_type")
    end

    # Verify that the user has provided the necessary job IDs
    isnothing(job_id) && error("Must pass GPU coupled run name to compare it.")

    # Construct artifacts directory
    output_dir = parsed_args["coupler_output_dir"]
    artifacts_dir = joinpath(output_dir, job_id, "artifacts")

    return (job_id, artifacts_dir)
end

"""
    append_table_data(table_data, setup_id, job_id, artifacts_dir)

Append data for a given setup to the table data.
"""
function append_table_data(table_data, setup_id, job_id, artifacts_dir)
    # Read in SYPD info
    sypd = open(joinpath(artifacts_dir, "sypd.txt"), "r") do sypd_file
        round(parse(Float64, read(sypd_file, String)), digits = 4)
    end

    # Create rows containing data for these runs
    new_table_data = [
        [setup_id "job ID:" job_id]
        ["" "SYPD:" sypd]
    ]
    return vcat(table_data, new_table_data)
end
