using Pkg
Pkg.activate(@__DIR__)

# First include and import the module
include("buildkite_analysis.jl")
using .BuildkiteAnalysis:
    set_offline_mode,
    analyze_buildkite_log,
    summarize_buildkite_errors,
    format_simulation_time

# Test configuration
pipeline = "climacoupler-longruns"
current_build_id = "915"
previous_build_id = "914"
job_id_step_names = [("amip_diagedmf_topo_integrated_land_gpu", "GPU AMIP + diag. EDMF + Earth topography + integrated land"),
("amip_edonly_topo_integrated_land_gpu", "GPU AMIP + ED only + Earth topography + integrated land")]

# Function to set up offline test environment with actual log file
function setup_offline_test()
    set_offline_mode(enabled = true, log_directory = "test_logs")
    mkpath("test_logs")

    # Path to the actual log file in Downloads
    downloads_log = joinpath(
        homedir(),
        "Downloads",
        "climaatmos-gpulongruns_build_577_computer-hydrostatic-balance.log",
    )

    if !isfile(downloads_log)
        error("Log file not found at: $downloads_log")
    end

    # Copy the log file to our test directory with the expected naming format
    cp(
        downloads_log,
        joinpath("test_logs", "build_577_job_computer-hydrostatic-balance.log"),
    )
end

# Function to clean up test files
function cleanup_test()
    rm("test_logs", recursive = true, force = true)
end

# Run the analysis with actual log file
function run_analysis(is_offline::Bool = true)
    is_offline && setup_offline_test()

    pipeline = "climacoupler-longruns"

    try
        summary = summarize_buildkite_errors(
            pipeline,
            current_build_id,
            previous_build_id,
            job_id_step_names,
        )
        println("\nBuildkite Analysis Summary:")
        println("=================================")
        println(summary)
    catch e
        println("Error during analysis: ", e)
        println(sprint(showerror, e, catch_backtrace()))
    finally
        cleanup_test()
    end
end

# Run the analysis
run_analysis(false)
