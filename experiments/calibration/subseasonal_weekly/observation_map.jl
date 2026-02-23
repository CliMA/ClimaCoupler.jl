using Dates: Week, Month, Year, Day, Millisecond
import JLD2
import ClimaAnalysis
import ClimaCoupler
import ClimaCalibrate.Checker: SequentialIndicesChecker

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

# Load calibration setup for NORMALIZE_VARIABLES setting
include(joinpath(@__DIR__, "calibration_setup.jl"))

# Load normalization stats if normalization is enabled
const NORM_STATS_PATH = joinpath(
    pkgdir(ClimaCoupler),
    "experiments/calibration/subseasonal_weekly/norm_stats.jld2",
)

# Global variable to cache normalization stats (loaded lazily)
const _NORM_STATS_CACHE = Ref{Union{Nothing, Dict{String, Tuple{Float64, Float64}}}}(nothing)

"""
    get_norm_stats()

Load and cache normalization statistics. Returns nothing if normalization is disabled
or stats file doesn't exist.
"""
function get_norm_stats()
    if !NORMALIZE_VARIABLES
        return nothing
    end
    
    if isnothing(_NORM_STATS_CACHE[])
        if isfile(NORM_STATS_PATH)
            _NORM_STATS_CACHE[] = JLD2.load_object(NORM_STATS_PATH)
            @info "Loaded normalization stats from $NORM_STATS_PATH"
        else
            @warn "Normalization enabled but norm_stats.jld2 not found at $NORM_STATS_PATH"
            @warn "Run generate_observations.jl first to create normalization stats"
            return nothing
        end
    end
    
    return _NORM_STATS_CACHE[]
end

"""
    normalize_model_var(var, norm_stats)

Normalize a model OutputVar using the same stats computed from observations.
"""
function normalize_model_var(var, norm_stats)
    short_name = ClimaAnalysis.short_name(var)
    if isnothing(norm_stats) || !haskey(norm_stats, short_name)
        return var
    end
    
    mean_val, std_val = norm_stats[short_name]
    normalized_data = (var.data .- mean_val) ./ std_val
    
    return ClimaAnalysis.OutputVar(
        var.attributes,
        var.dims,
        var.dim_attributes,
        normalized_data,
    )
end

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


# For single-day/weekly calibration, we just need to set units and window by time
# observations are already resampled to model grid in generate_observations.jl
function preprocess_var(var, sample_date_range)
    # Set units to match what's expected
    var = set_units(var, var_units[ClimaAnalysis.short_name(var)])
    # Window to the date range
    var = ClimaAnalysis.window(var, "time"; left = sample_date_range[1], right = sample_date_range[2])
    return var
end

# Override process_member_data! to apply normalization to model output
function process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx, iteration)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges[iteration + 1]
    @info "Short names: $short_names"

    # Load normalization stats (cached after first call)
    norm_stats = get_norm_stats()
    if !isnothing(norm_stats)
        @info "Applying normalization to model output"
    end

    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in short_names
        var = get_var(short_name, simdir)
        var = preprocess_var(var, sample_date_ranges)
        
        # Apply normalization if enabled
        if !isnothing(norm_stats)
            var = normalize_model_var(var, norm_stats)
        end

        EnsembleBuilder.fill_g_ens_col!(
            g_ens_builder,
            col_idx,
            var;
            checkers = (SequentialIndicesChecker(),),
            verbose = true,
        )
    end

    return nothing
end

# Get job_id from config file name (e.g., "wxquest_diagedmf_weekly_calibration.yml" -> "wxquest_diagedmf_weekly_calibration")
function get_job_id()
    config_file = CALIBRATE_CONFIG.config_file
    return replace(basename(config_file), ".yml" => "")
end

# Override observation_map to use correct job_id path
function ClimaCalibrate.observation_map(iteration)
    output_dir = CALIBRATE_CONFIG.output_dir
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))

    g_ens_builder = EnsembleBuilder.GEnsembleBuilder(ekp)
    job_id = get_job_id()

    for m in 1:EKP.get_N_ens(ekp)
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, job_id, "output_active")
        @info "Processing member $m: $simdir_path"
        try
            process_member_data!(g_ens_builder, simdir_path, m, iteration)
        catch e
            @error "Ensemble member $m failed" exception = (e, catch_backtrace())
        end
    end

    g_ensemble = EnsembleBuilder.get_g_ensemble(g_ens_builder)

    # Too many NaNs - abort
    if sum(isnan.(g_ensemble)) > 0.5 * length(g_ensemble)
        error("Too many NaNs")
    end

    return g_ensemble
end

# Override analyze_iteration (placeholder for now)
function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    @info "Iteration $iteration analysis (placeholder)"
    return nothing
end
