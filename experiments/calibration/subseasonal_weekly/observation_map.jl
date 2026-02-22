using Dates: Week, Month, Year, Day, Millisecond
import JLD2

# Override JLD2's default_iotype to use IOStream instead of MmapIO
# This avoids Bus errors from memory-mapped files on Lustre filesystem
JLD2.default_iotype() = IOStream

# Reuse observation_map, process_member_data!, etc. from subseasonal
include(
    joinpath(
        @__DIR__,
        "..",
        "subseasonal",
        "observation_map.jl",
    ),
)

# Override largest_period with weekly-specific thresholds
function largest_period(sample_date_range)
    span = maximum(sample_date_range) - minimum(sample_date_range)
    span = Millisecond(span)
    # For weekly data (6 days span = 7 day period including endpoints)
    span.value == 0 && return Week(1)
    day_in_ms = 8.64e7
    period =
        span.value >= day_in_ms * 365 ? Year(1) :
        span.value >= day_in_ms * 28 ? Month(1) :
        span.value >= day_in_ms * 6 ? Week(1) : Day(1)
    return period
end

# Override analyze_iteration (placeholder for now)
function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    @info "Iteration $iteration analysis (placeholder)"
    return nothing
end
