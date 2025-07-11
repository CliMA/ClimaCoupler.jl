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
current_build_id = "577"
previous_build_id = "576"
job_ids = [("computer-hydrostatic-balance", "GPU Long Run Test")]

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
function run_analysis()
    setup_offline_test()

    try
        summary = summarize_buildkite_errors(
            current_build_id,
            previous_build_id,
            job_ids,
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
run_analysis()
