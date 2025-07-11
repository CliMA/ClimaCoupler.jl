module BuildkiteAnalysis

import HTTP
import JSON
using Dates
using Logging

# Export public functions
export set_offline_mode,
    analyze_buildkite_log,
    summarize_buildkite_errors,
    BuildkiteError,
    format_simulation_time

"""
    BuildkiteError
Structure containing error information from a Buildkite job
"""
struct BuildkiteError
    error_code::Union{Int, Nothing}
    last_sim_time::Union{Float64, Nothing}
    job_id::String
    build_id::String
    step_name::String
    timestamp::DateTime
    clima_job_id::String  # ClimaAtmos job identifier
    sypd::Union{Float64, Nothing}  # Simulated Years Per Day
    build_status::String  # Success/Failure/Cancelled status
end

# Add offline mode configuration
mutable struct Config
    offline_mode::Bool
    log_directory::String
end

# Global configuration
const GLOBAL_CONFIG = Config(false, "buildkite_logs")

"""
    set_offline_mode(;enabled::Bool = true, log_directory::String = "buildkite_logs")
Configure the module to run in offline mode using local log files.
"""
function set_offline_mode(;
    enabled::Bool = true,
    log_directory::String = "buildkite_logs",
)
    GLOBAL_CONFIG.offline_mode = enabled
    GLOBAL_CONFIG.log_directory = log_directory
end

"""
    save_log_locally(build_id::String, job_id::String, content::String)
Save a Buildkite log to a local file for offline use.
"""
function save_log_locally(build_id::String, job_id::String, content::String)
    # Create log directory if it doesn't exist
    mkpath(GLOBAL_CONFIG.log_directory)

    # Create a filename based on build and job IDs
    filename = joinpath(
        GLOBAL_CONFIG.log_directory,
        "build_$(build_id)_job_$(job_id).log",
    )

    # Save the content
    write(filename, content)
    @info "Saved log to $filename"
end

"""
    read_log_locally(build_id::String, job_id::String)
Read a Buildkite log from a local file.
"""
function read_log_locally(build_id::String, job_id::String)
    filename = joinpath(
        GLOBAL_CONFIG.log_directory,
        "build_$(build_id)_job_$(job_id).log",
    )
    if !isfile(filename)
        error("Log file not found: $filename")
    end
    return read(filename, String)
end

"""
    extract_sim_time(log_text::String)
Extract the simulation time from ClimaAtmos progress reports.
Returns the time in seconds if found, nothing otherwise.
"""
function extract_sim_time(log_text::String)
    # Look for the progress report format
    pattern = r"simulation_time = \"(\d+) weeks?, (\d+) days?\""

    # Find all matches and take the last one (most recent)
    matches = collect(eachmatch(pattern, log_text))
    if !isempty(matches)
        last_match = last(matches)
        weeks = parse(Int, last_match[1])
        days = parse(Int, last_match[2])

        # Convert to seconds
        total_days = weeks * 7 + days
        return total_days * 86400.0  # days to seconds
    end

    return nothing
end

"""
    extract_error_code(log_text::String)
Extract error code from log text.
Returns nothing if no error code is found.
"""
function extract_error_code(log_text::String)
    # Common patterns in Buildkite logs
    patterns = [
        r"Error: (\d+)",
        r"exit code (\d+)",
        r"Exit status: (\d+)",
        r"Process exited with status (\d+)",
        r"ERROR: Process completed with exit code (\d+)",
        r"Test Failed at .* status (\d+)",
        r"\[error\] Exit code: (\d+)",
    ]

    for pattern in patterns
        match_result = match(pattern, log_text)
        if !isnothing(match_result)
            return parse(Int, match_result[1])
        end
    end

    # If we find specific error patterns without codes, assign custom codes
    error_patterns = [
        (r"ERROR: LoadError:", 1),
        (r"signal \((\d+)\)", -1),  # Capture the signal number as negative
        (r"ERROR: SystemError", 2),
        (r"ERROR: ArgumentError", 3),
        (r"Test Failed at", 1),
        (r"ERROR: AssertionError", 4),
    ]

    for (pattern, code) in error_patterns
        if match(pattern, log_text) !== nothing
            return code
        end
    end

    # Check for timeout or cancellation
    if contains(log_text, "Cancelling") || contains(log_text, "Canceled")
        return -15  # Common signal for termination
    end

    return nothing
end

"""
    fetch_buildkite_log(build_id::String, job_id::String)
Fetch log content from Buildkite API or local file system depending on mode.
"""
function fetch_buildkite_log(build_id::String, job_id::String)
    if GLOBAL_CONFIG.offline_mode
        log_file = joinpath(
            GLOBAL_CONFIG.log_directory,
            "build_$(build_id)_job_$(job_id).log",
        )
        @info "Reading log for specific job" build_id job_id log_file
        return read_log_locally(build_id, job_id)
    else
        token = get(ENV, "BUILDKITE_TOKEN", nothing)
        if isnothing(token)
            error(
                "BUILDKITE_TOKEN environment variable must be set when not in offline mode",
            )
        end

        # Note: This API endpoint is specifically for a single job's logs
        api_url = "https://api.buildkite.com/v2/organizations/clima/pipelines/climaatmos-gpulongruns/builds/$build_id/jobs/$job_id/log"
        @info "Fetching log for specific job" build_id job_id api_url

        headers = Dict("Authorization" => "Bearer $token")
        response = HTTP.get(api_url, headers)
        content = String(response.body)

        # Save the job-specific log locally for future offline use
        try
            save_log_locally(build_id, job_id, content)
        catch e
            @warn "Failed to save log locally" build_id job_id exception=e
        end

        return content
    end
end

"""
    extract_job_id(log_text::String)
Extract the ClimaAtmos job identifier from the log.
Returns the job_id string if found, "unknown" otherwise.
"""
function extract_job_id(log_text::String)
    # Look for the job_id pattern in the log
    pattern = r"job_id = \"([^\"]+)\""
    match_result = match(pattern, log_text)

    if !isnothing(match_result)
        return match_result[1]
    end

    return "unknown"
end

"""
    extract_sypd(log_text::String)
Extract the estimated SYPD (Simulated Years Per Day) from the log.
Returns the SYPD as a float if found, nothing otherwise.
"""
function extract_sypd(log_text::String)
    # Look for the SYPD pattern in the log
    pattern = r"estimated_sypd = \"([0-9.]+)\""
    matches = collect(eachmatch(pattern, log_text))

    if !isempty(matches)
        # Take the last reported SYPD value
        return parse(Float64, last(matches)[1])
    end

    return nothing
end

"""
    extract_build_status(log_text::String)
Determine if the build was successful, failed, or cancelled based on log patterns.
Returns "Success", "Failed", or "Cancelled".
"""
function extract_build_status(log_text::String)
    # First verify we're looking at a single job's log
    if !occursin(r"job_id = \"[^\"]+\"", log_text)
        @warn "Log content might not be for a specific job"
    end

    # First check for cancellation patterns (highest priority)
    cancel_patterns = [
        r"Cancelling",
        r"Canceled",
        r"Build was canceled",
        r"Terminated",
        r"signal \(15\)",  # SIGTERM
    ]

    # Check for critical failure patterns
    critical_failure_patterns = [
        r"Test Failed at",
        r"Process exited with status [1-9]",
        r"Build Failed",
        r"Error running command",
        r"signal \(9\)",  # SIGKILL
        r"ERROR: LoadError:",
    ]

    # Check for success patterns
    success_patterns = [
        r"All tests pass",
        r"Process exited with status 0",
        r"Build Finished",
        r"Simulation completed successfully",
        r"percent_complete = \"100%\"",
        r"percent_complete = \"[9][5-9]\.?\d*%\"",  # 95-99.x%
        r"estimated_finish_date",  # Progress report indicates running simulation
    ]

    # First check for cancellation
    for pattern in cancel_patterns
        if occursin(pattern, log_text)
            return "Cancelled"
        end
    end

    # Then check for critical failures
    for pattern in critical_failure_patterns
        if occursin(pattern, log_text)
            return "Failed"
        end
    end

    # Check for success indicators
    for pattern in success_patterns
        if occursin(pattern, log_text)
            # Additional verification: if we find a success pattern,
            # make sure there's no critical failure after it
            match_pos = findfirst(pattern, log_text)
            if !isnothing(match_pos)
                # Only check the log after this success pattern
                remaining_log = log_text[match_pos.stop:end]
                has_failure = any(
                    p -> occursin(p, remaining_log),
                    critical_failure_patterns,
                )
                if !has_failure
                    return "Success"
                end
            end
        end
    end

    # For ClimaAtmos, if we see progress reports with high completion percentage
    # and no critical failures, consider it a success
    if occursin(r"percent_complete = \"[0-9]+\.?[0-9]*%\"", log_text) &&
       !any(p -> occursin(p, log_text), critical_failure_patterns)
        return "Success"
    end

    # If we have simulation progress and no critical failures, consider it running successfully
    if occursin(r"simulation_time = \".*\"", log_text) &&
       !any(p -> occursin(p, log_text), critical_failure_patterns)
        return "Success"
    end

    # Default to Failed only if we can't determine success
    return "Failed"
end

"""
    analyze_buildkite_log(build_id::String, job_id::String, step_name::String)
Analyze a single Buildkite job log and return structured error information.
"""
function analyze_buildkite_log(
    build_id::String,
    job_id::String,
    step_name::String,
)
    @info "Analyzing specific job log" build_id job_id step_name

    log_content = fetch_buildkite_log(build_id, job_id)

    # Verify we have the correct job's log by checking for job_id in the content
    if !isnothing(match(r"job_id = \"([^\"]+)\"", log_content))
        extracted_job_id = extract_job_id(log_content)
        @info "Found job identifier in log" extracted_job_id
    end

    error_code = extract_error_code(log_content)
    sim_time = extract_sim_time(log_content)
    clima_job_id = extract_job_id(log_content)
    sypd = extract_sypd(log_content)
    build_status = extract_build_status(log_content)

    @info "Analysis results for job" build_id job_id build_status error_code

    return BuildkiteError(
        error_code,
        sim_time,
        job_id,
        build_id,
        step_name,
        now(),
        clima_job_id,
        sypd,
        build_status,
    )
end

"""
    format_simulation_time(seconds::Float64)
Convert simulation time to weeks and days format as used in ClimaAtmos.
"""
function format_simulation_time(seconds::Float64)
    total_days = seconds / 86400  # Convert seconds to days
    weeks = floor(Int, total_days / 7)
    remaining_days = round(Int, total_days % 7)

    # Format string components
    week_str = weeks == 1 ? "week" : "weeks"
    day_str = remaining_days == 1 ? "day" : "days"

    if weeks > 0
        return "$(weeks) $(week_str), $(remaining_days) $(day_str)"
    else
        return "$(remaining_days) $(day_str)"
    end
end

"""
    generate_build_report(current::BuildkiteError, previous::Union{BuildkiteError, Nothing}=nothing)
Generate a report comparing two builds, or just current build if previous is nothing.
"""
function generate_build_report(
    current::BuildkiteError,
    previous::Union{BuildkiteError, Nothing} = nothing,
)
    status_emoji = Dict("Success" => "✅", "Failed" => "❌", "Cancelled" => "⚠️")

    report = String[]

    # Build Information
    push!(report, "Build Information:")
    push!(report, "─" * repeat("─", 30))
    push!(report, "ClimaAtmos Job ID: $(current.clima_job_id)")
    push!(
        report,
        "Current Build: $(current.build_id) $(get(status_emoji, current.build_status, "❓")) ($(current.build_status))",
    )
    push!(
        report,
        "Current Build Date: $(Dates.format(current.timestamp, "yyyy-mm-dd HH:MM:SS"))",
    )
    push!(
        report,
        isnothing(previous) ? "Previous Build: Not Available" :
        "Previous Build: $(previous.build_id) $(get(status_emoji, previous.build_status, "❓")) ($(previous.build_status))",
    )
    push!(
        report,
        isnothing(previous) ? "Previous Build Date: Not Available" :
        "Previous Build Date: $(Dates.format(previous.timestamp, "yyyy-mm-dd HH:MM:SS"))",
    )
    push!(report, "")

    # Simulation Progress
    push!(report, "Simulation Progress:")
    push!(report, "─" * repeat("─", 30))
    push!(
        report,
        "Current Run: $(isnothing(current.last_sim_time) ? "No simulation time found" : format_simulation_time(current.last_sim_time))",
    )

    if !isnothing(previous)
        push!(
            report,
            "Previous Run: $(isnothing(previous.last_sim_time) ? "No simulation time found" : format_simulation_time(previous.last_sim_time))",
        )
        if !isnothing(current.last_sim_time) &&
           !isnothing(previous.last_sim_time)
            time_diff = current.last_sim_time - previous.last_sim_time
            diff_days = abs(time_diff) / 86400
            diff_str =
                diff_days >= 7 ?
                "$(floor(Int, diff_days/7)) weeks, $(round(Int, diff_days%7)) days" :
                "$(round(Int, diff_days)) days"
            push!(
                report,
                "Progress Difference: $(diff_str) ($(time_diff < 0 ? "behind" : "ahead"))",
            )
        end
    else
        push!(report, "Previous Run: Not Available")
    end
    push!(report, "")

    # Performance Metrics
    push!(report, "Performance Metrics:")
    push!(report, "─" * repeat("─", 30))
    push!(
        report,
        "Current Run SYPD: $(isnothing(current.sypd) ? "Not available" : round(current.sypd, digits=3))",
    )

    if !isnothing(previous)
        push!(
            report,
            "Previous Run SYPD: $(isnothing(previous.sypd) ? "Not available" : round(previous.sypd, digits=3))",
        )
        if !isnothing(current.sypd) && !isnothing(previous.sypd)
            sypd_diff = current.sypd - previous.sypd
            sypd_change_pct = (sypd_diff / previous.sypd) * 100
            push!(
                report,
                "SYPD Change: $(round(sypd_diff, digits=3)) ($(round(sypd_change_pct, digits=1))%)",
            )
        end
    else
        push!(report, "Previous Run SYPD: Not Available")
    end
    push!(report, "")

    # Error Information
    push!(report, "Error Information:")
    push!(report, "─" * repeat("─", 30))
    push!(
        report,
        "Current Run Error Code: $(isnothing(current.error_code) ? "None" : current.error_code)",
    )
    if !isnothing(previous)
        push!(
            report,
            "Previous Run Error Code: $(isnothing(previous.error_code) ? "None" : previous.error_code)",
        )
        if !isnothing(current.error_code) &&
           !isnothing(previous.error_code) &&
           current.error_code != previous.error_code
            push!(
                report,
                "Error Code Changed: $(previous.error_code) → $(current.error_code)",
            )
        end
    else
        push!(report, "Previous Run: Not Available")
    end

    return join(report, "\n")
end

"""
    summarize_buildkite_errors(current_build::String, previous_build::String, job_ids::Vector{Tuple{String,String}})
Generate summary reports for multiple jobs, comparing current and previous builds.
"""
function summarize_buildkite_errors(
    current_build::String,
    previous_build::String,
    job_ids::Vector{Tuple{String, String}},
)
    # Temporarily disable info logging
    current_logger = Logging.global_logger()
    Logging.global_logger(Logging.SimpleLogger(stderr, Logging.Warn))

    try
        summaries = String[]
        for (job_id, step_name) in job_ids
            try
                current_error =
                    analyze_buildkite_log(current_build, job_id, step_name)
                previous_error = try
                    analyze_buildkite_log(previous_build, job_id, step_name)
                catch e
                    isa(e, SystemError) ||
                        contains(string(e), "Log file not found") ?
                    nothing : rethrow(e)
                end
                push!(
                    summaries,
                    generate_build_report(current_error, previous_error),
                )
            catch e
                push!(
                    summaries,
                    "Failed to analyze $step_name: $(sprint(showerror, e))",
                )
            end
        end
        return join(summaries, "\n\n")
    finally
        # Restore original logger
        Logging.global_logger(current_logger)
    end
end

end # module 
