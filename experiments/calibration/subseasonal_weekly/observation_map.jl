using Dates: Week, Month, Year, Day, Millisecond
import JLD2
import ClimaAnalysis
import ClimaCoupler
import ClimaCalibrate.Checker: SequentialIndicesChecker
import CairoMakie
import GeoMakie

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


# For monthly calibration, window model output to the calibration period
# The spinup is already handled by model_interface.jl (model starts earlier)
# sample_date_range = (calib_start, calib_end) defines the comparison period
# Then TIME-AVERAGE to a single value for comparison to CERES monthly mean
function preprocess_var(var, sample_date_range)
    # Set units to match what's expected
    var = set_units(var, get_var_units(ClimaAnalysis.short_name(var)))
    
    # Window to calibration period defined by sample_date_range
    calib_start, calib_end = sample_date_range
    
    @info "Windowing model output: $calib_start to $calib_end"
    var = ClimaAnalysis.window(var, "time"; left = calib_start, right = calib_end)
    
    # Time-average to single value for comparison to CERES monthly mean
    n_times = length(ClimaAnalysis.times(var))
    if n_times > 1
        @info "Time-averaging $n_times daily outputs to single mean"
        var = time_average_with_date(var, calib_start)
    end
    
    return var
end

"""
    time_average_with_date(var, date)

Average over time dimension but preserve a single time point at `date`.
This is needed because ClimaCalibrate expects OutputVars with time dimensions
to match observation dates.
"""
function time_average_with_date(var, date)
    # Get the time-averaged data (2D: lon x lat, time dimension removed)
    avg_var = ClimaAnalysis.average_time(var)
    
    # Set start_date to the target date and time_val to 0
    # This ensures ClimaAnalysis.dates() returns [date], matching the observation
    start_date_str = Dates.format(date, "yyyy-mm-dd")
    new_attribs = Dict{String, Any}(k => v for (k, v) in avg_var.attributes)
    new_attribs["start_date"] = start_date_str
    # Create new dims with time added back
    new_dims = copy(avg_var.dims)
    new_dims["time"] = [0.0]
    new_dim_attribs = copy(avg_var.dim_attributes)
    new_dim_attribs["time"] = Dict{String, Any}("units" => "s")
    # Reshape data to add time dimension (lon x lat -> lon x lat x 1)
    new_data = reshape(avg_var.data, size(avg_var.data)..., 1)
    
    return ClimaAnalysis.OutputVar(new_attribs, new_dims, new_dim_attribs, new_data)
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
            # Fill failed member column with NaN so EKP can handle the failure
            EnsembleBuilder.fill_g_ens_col!(g_ens_builder, m, NaN)
        end
    end

    g_ens = EnsembleBuilder.get_g_ensemble(g_ens_builder)
    # Too many NaNs - abort (90% threshold like subseasonal)
    if count(isnan, g_ens) > 0.9 * length(g_ens)
        error("Too many NaNs")
    end
    return EnsembleBuilder.is_complete(g_ens_builder) ? g_ens :
           error("G ensemble matrix is not completed")
end

# Extend bias_plot_extrema (defined in subseasonal/observation_map.jl) with weekly variables
bias_plot_extrema["ta_850hPa"] = (-2, 2)
bias_plot_extrema["ta_500hPa"] = (-2, 2)
bias_plot_extrema["ta_200hPa"] = (-2, 2)
bias_plot_extrema["hur_850hPa"] = (-2, 2)
bias_plot_extrema["hur_500hPa"] = (-2, 2)
bias_plot_extrema["hur_200hPa"] = (-2, 2)

# TODO: Unify this with `plot_bias` in the subseasonal observation_map.jl
"""
    plot_bias_weekly(ekp, simdir, iteration; output_dir)

Plot bias maps comparing simulation output to ERA5 observations for all variables in
`CALIBRATE_CONFIG.short_names`. ERA5 vars are reconstructed from the EKP observation
object and denormalized to physical units when normalization is enabled.
"""
function plot_bias_weekly(ekp, simdir, iteration; output_dir = simdir.simulation_path)
    sample_date_range = CALIBRATE_CONFIG.sample_date_ranges[iteration + 1]
    calib_start, _ = sample_date_range

    # Reconstruct ERA5 OutputVars from the EKP observation object
    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs = ClimaCalibrate.ObservationRecipe.get_observations_for_nth_iteration(
        obs_series,
        iteration + 1,
    )

    # ERA5 vars are kept as-is (normalized), matching what enters the loss function
    era5_vars = mapreduce(ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, minibatch_obs)

    # Normalize sim vars identically to process_member_data! so the bias reflects
    # exactly the residual seen by the loss function
    norm_stats = get_norm_stats()
    sim_vars = map(CALIBRATE_CONFIG.short_names) do short_name
        var = preprocess_var(get_var(short_name, simdir), sample_date_range)
        normalize_model_var(var, norm_stats)
    end

    # Match sim_vars with era5_vars by short_name
    var_pairs = []
    for sim_var in sim_vars
        sn = ClimaAnalysis.short_name(sim_var)
        era5_idx = findfirst(v -> ClimaAnalysis.short_name(v) == sn, era5_vars)
        if !isnothing(era5_idx)
            push!(var_pairs, (sim_var, era5_vars[era5_idx]))
        else
            @warn "No ERA5 data found for $sn — skipping bias plot"
        end
    end

    if isempty(var_pairs)
        @warn "No matching variable pairs found for bias plotting"
        return nothing
    end

    fig = GeoMakie.Figure(size = (1500, 500 * length(var_pairs)))
    for (i, (sim_var, era5_var)) in enumerate(var_pairs)
        sn = ClimaAnalysis.short_name(sim_var)
        sim_var_t = ClimaAnalysis.slice(sim_var, time = calib_start)
        era5_var_t = ClimaAnalysis.slice(era5_var, time = calib_start)
        cmap_extrema = get(bias_plot_extrema, sn, extrema(sim_var_t.data))
        ClimaAnalysis.Visualize.plot_bias_on_globe!(
            fig[i, 1],
            sim_var_t,
            era5_var_t;
            cmap_extrema,
        )
    end

    GeoMakie.save(joinpath(output_dir, "bias_sample_dates.png"), fig)
    return nothing
end

function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
    plot_constrained_params_and_errors(output_dir, ekp, prior)

    job_id = get_job_id()
    member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1)
    simdir_path = joinpath(member_path, job_id, "output_active")
    try
        simdir = ClimaAnalysis.SimDir(simdir_path)
        plot_bias_weekly(ekp, simdir, iteration; output_dir = plot_output_path)
    catch e
        @error "Bias plotting failed" exception = (e, catch_backtrace())
    end

    @info "Ensemble spread: $(scalar_spread(ekp))"
    return nothing
end
