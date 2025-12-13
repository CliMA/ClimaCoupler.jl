#=
Our goal here is to output a table displaying some results from benchmark runs
in the coupler. We mainly want to compare SYPD between coupled and atmos-only runs.

The table should look something like this:
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

import PrettyTables
import ClimaCoupler
import ClimaCoupler.SimOutput

# Read in command line arguments
parsed_args = SimOutput.get_benchmark_args()
output_dir = parsed_args["coupler_output_dir"]

# Access buildkite pipeline ID (from `BUILDKITE_GITHUB_DEPLOYMENT_ID` variable)
build_id = parsed_args["build_id"]
if !isnothing(build_id)
    build_id_str = "Build ID: $build_id"
else
    build_id_str = "Build ID: N/A"
end

# Set up info for PrettyTables.jl
headers = [build_id_str, "Horiz. res.: 30 elems", "GPU Run [2 A100s]"]
data = [
    ["" "Vert. res.: 63 levels" ""]
    ["" "dt: 120secs" ""]
]

# Set up a list of tuples containing the run name and description
run_names = [
    ("coupled_progedmf_coarse", "Coupled with progedmf + 1M (16 helem)"),
    ("coupled_progedmf_fine", "Coupled with progedmf + 1M (30 helem)"),
    ("coupled_io", "Coupled with diag. EDMF + IO"),
    ("coupled", "Coupled with diag. EDMF"),
    ("atmos_diagedmf", "Atmos with diag. EDMF"),
    ("atmos", "Atmos without diag. EDMF"),
]

# For each run, get the run info and append it to the table
for (run_name, description) in run_names
    run_info = SimOutput.get_run_info(parsed_args, run_name)
    data = SimOutput.append_table_data(data, description, run_info...)
end

# Output table to file, note that this must match the slack upload command in the pipeline.yml file
table_output_dir = joinpath(output_dir, "compare_amip_climaatmos_amip_diagedmf")
mkpath(table_output_dir)
open(joinpath(table_output_dir, "table.txt"), "w") do f
    # Output the table, including lines before and after the header
    PrettyTables.pretty_table(
        f,
        data,
        header = headers,
        hlines = [0, 3, 5, 7, 9, 11, 13, 15],
    )
end
