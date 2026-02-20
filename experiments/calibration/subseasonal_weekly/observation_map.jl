using Statistics
import JLD2
import Dates
using Dates: Week, Month, Year, Day, Millisecond
using ClimaAnalysis
import ClimaCalibrate
import ClimaCoupler
import ClimaCalibrate: EnsembleBuilder
import EnsembleKalmanProcesses as EKP

# Override JLD2's default_iotype to use IOStream instead of MmapIO
# This avoids Bus errors from memory-mapped files on Lustre filesystem
JLD2.default_iotype() = IOStream

include(joinpath(@__DIR__, "observation_utils.jl"))

"""
    retry_on_error(f; max_retries=5, delay=2.0, backoff=2.0)

Retry a function `f` on any error, with exponential backoff.
Useful for handling transient filesystem errors on Lustre (NetCDF/HDF5 race conditions).
"""
function retry_on_error(f; max_retries = 5, delay = 2.0, backoff = 2.0)
    last_error = nothing
    for attempt in 1:max_retries
        try
            return f()
        catch e
            last_error = e
            if attempt < max_retries
                wait_time = delay * backoff^(attempt - 1)
                @warn "Attempt $attempt failed, retrying in $(wait_time)s..." exception =
                    (e, catch_backtrace())
                sleep(wait_time)
            end
        end
    end
    @error "All $max_retries attempts failed"
    throw(last_error)
end

"""
    ClimaCalibrate.observation_map(iteration)

Return G ensemble for an `iteration` using EnsembleBuilder.

G ensemble represents the concatenated forward model evaluations from all
ensemble members, arranged horizontally. Each individual forward model
evaluation corresponds to preprocessed, flattened simulation data from a single
ensemble member that has been matched to the corresponding observational data.

Uses EnsembleBuilder to ensure G vector format matches observation vector format
(created by ObservationRecipe in generate_observations.jl).
"""
function ClimaCalibrate.observation_map(iteration)
    output_dir = CALIBRATE_CONFIG.output_dir
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))

    g_ens_builder = EnsembleBuilder.GEnsembleBuilder(ekp)

    # Small delay to ensure all files are flushed to Lustre before reading
    sleep(2.0)

    for m in 1:EKP.get_N_ens(ekp)
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        subdirs = filter(isdir, [joinpath(member_path, f) for f in readdir(member_path)])
        if isempty(subdirs)
            error("No output directory found in $member_path")
        end
        job_dir = first(subdirs)
        simdir_path = joinpath(job_dir, "output_active", "clima_atmos")
        @info "Processing member $m: $simdir_path"
        try
            process_member_data!(g_ens_builder, simdir_path, m, iteration)
        catch e
            @error "Ensemble member $m failed" exception = (e, catch_backtrace())
            EnsembleBuilder.fill_g_ens_col!(g_ens_builder, m, NaN)
        end
    end

    g_ens = EnsembleBuilder.get_g_ensemble(g_ens_builder)
    if count(isnan, g_ens) > 0.9 * length(g_ens)
        error("Too many NaNs in G ensemble")
    end
    return EnsembleBuilder.is_complete(g_ens_builder) ? g_ens :
           error("G ensemble matrix is not completed")
end

"""
    process_member_data!(g_ens_builder, simdir_path, col_idx, iteration)

Process a single ensemble member's simulation output and fill the corresponding
column of the G ensemble matrix using EnsembleBuilder.

This mirrors the subseasonal pipeline pattern: get vars, preprocess, and pass
to `fill_g_ens_col!` which handles flattening/alignment to match the observation
vector format (created by ObservationRecipe in generate_observations.jl).
"""
function process_member_data!(g_ens_builder, simdir_path, col_idx, iteration)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    sample_date_range = first(CALIBRATE_CONFIG.sample_date_ranges)

    # Wrap SimDir in retry logic for transient Lustre/NetCDF errors
    simdir = retry_on_error() do
        ClimaAnalysis.SimDir(simdir_path)
    end

    for short_name in short_names
        var = retry_on_error() do
            get_var(short_name, simdir)
        end
        var = preprocess_var(var, sample_date_range)

        EnsembleBuilder.fill_g_ens_col!(
            g_ens_builder,
            col_idx,
            var;
            checkers = (EnsembleBuilder.SequentialIndicesChecker(),),
            verbose = true,
        )
    end

    return nothing
end

# get_var and preprocess_var are imported from subseasonal/observation_utils.jl (shared).
# largest_period is overridden here with weekly-specific thresholds.

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

"""
    ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)

Placeholder for iteration analysis (bias plots, etc.)
"""
function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    @info "Iteration $iteration analysis (placeholder)"
    return nothing
end
