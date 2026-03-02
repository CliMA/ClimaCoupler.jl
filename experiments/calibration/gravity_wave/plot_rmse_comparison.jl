#=
plot_rmse_comparison.jl

Compare initial vs calibrated ensemble bias against ERA5 observations.

Usage:
    julia --project=experiments/ClimaEarth plot_rmse_comparison.jl [output_dir]

If output_dir is not provided, defaults to "output/gw_calibration"
=#

using ClimaAnalysis
using CairoMakie
using GeoMakie
using JLD2
using Statistics
using Printf

# Override JLD2's default_iotype to avoid Lustre issues
JLD2.default_iotype() = IOStream

# ==========================================================================
# CONFIGURATION
# ==========================================================================
output_dir = length(ARGS) >= 1 ? ARGS[1] : "output/gw_calibration"
output_dir = abspath(output_dir)

# Create plots subdirectory
plots_dir = joinpath(output_dir, "plots")
mkpath(plots_dir)

exp_suffix = basename(output_dir)

@info "RMSE Comparison Plot"
@info "Output directory: $output_dir"

# Variables to plot
short_names = ["tas", "mslp", "pr"]

# Units for display
var_units = Dict(
    "tas" => "K",
    "mslp" => "Pa",
    "pr" => "kg m⁻² s⁻¹",
)

# Colorbar ranges for bias plots
bias_ranges = Dict(
    "tas" => (-6, 6),
    "mslp" => (-1000, 1000),
    "pr" => (-1e-4, 1e-4),
)

# Colorbar ranges for improvement plots
improvement_ranges = Dict(
    "tas" => (-3, 3),
    "mslp" => (-500, 500),
    "pr" => (-5e-5, 5e-5),
)

# ==========================================================================
# DATA LOADING
# ==========================================================================

"""
    load_era5_var(short_name)

Load preprocessed ERA5 variable from gw_wq folder.
Returns the first matching variable (by short_name key pattern).
"""
function load_era5_var(short_name)
    # Try loading from gw_wq preprocessed vars
    preprocessed_path = joinpath(
        dirname(dirname(@__DIR__)),
        "calibration/gw_wq/preprocessed_vars.jld2"
    )

    if !isfile(preprocessed_path)
        error("Preprocessed vars not found at $preprocessed_path")
    end

    preprocessed_vars = JLD2.load_object(preprocessed_path)

    # Find the key matching this short_name
    for (key, var) in preprocessed_vars
        if key[1] == short_name
            @info "  Loaded ERA5 $short_name from $(key)"
            return var
        end
    end

    error("No ERA5 var found for $short_name")
end

"""
    load_ensemble_mean(output_dir, iteration, short_name)

Load simulation outputs for all members at given iteration and compute ensemble mean.
"""
function load_ensemble_mean(output_dir, iteration, short_name)
    iter_str = @sprintf("iteration_%03d", iteration)
    iter_path = joinpath(output_dir, iter_str)

    # Find all member directories
    member_dirs = filter(readdir(iter_path)) do name
        startswith(name, "member_") && isdir(joinpath(iter_path, name))
    end

    if isempty(member_dirs)
        error("No member directories found in $iter_path")
    end

    @info "Loading $short_name from $(length(member_dirs)) members at iteration $iteration"

    member_vars = []
    for member_dir in sort(member_dirs)
        member_path = joinpath(iter_path, member_dir)

        # Find the job output folder
        subdirs = filter(isdir, [joinpath(member_path, f) for f in readdir(member_path)])
        if isempty(subdirs)
            @warn "No output directory in $member_path, skipping"
            continue
        end
        job_dir = first(subdirs)
        simdir_path = joinpath(job_dir, "output_active", "clima_atmos")

        if !isdir(simdir_path)
            @warn "No clima_atmos directory in $job_dir, skipping"
            continue
        end

        try
            simdir = ClimaAnalysis.SimDir(simdir_path)
            var = ClimaAnalysis.get(simdir; short_name)

            # Average over time if multiple time steps
            if "time" in keys(var.dims)
                n_times = length(ClimaAnalysis.times(var))
                if n_times > 1
                    var = ClimaAnalysis.average_time(var)
                else
                    var = ClimaAnalysis.slice(var; time = first(ClimaAnalysis.times(var)))
                end
            end

            push!(member_vars, var)
        catch e
            @warn "Failed to load $short_name from $simdir_path" exception=e
        end
    end

    if isempty(member_vars)
        error("No valid member data for $short_name at iteration $iteration")
    end

    # Compute ensemble mean
    ensemble_mean = member_vars[1]
    for i in 2:length(member_vars)
        ensemble_mean = ensemble_mean + member_vars[i]
    end
    ensemble_mean = ensemble_mean / length(member_vars)

    @info "  Computed ensemble mean from $(length(member_vars)) members"
    return ensemble_mean
end

"""
    compute_global_rmse(sim_var, obs_var)

Compute global RMSE (weighted by latitude).
"""
function compute_global_rmse(sim_var, obs_var)
    diff = sim_var - obs_var
    sq_diff = ClimaAnalysis.replace(x -> x^2, diff)
    mean_sq = ClimaAnalysis.average_lonlat(sq_diff; weighted=true)
    return sqrt(mean_sq.data[1])
end

"""
    compute_global_bias(sim_var, obs_var)

Compute global mean bias (weighted by latitude).
"""
function compute_global_bias(sim_var, obs_var)
    diff = sim_var - obs_var
    mean_diff = ClimaAnalysis.average_lonlat(diff; weighted=true)
    return mean_diff.data[1]
end

# ==========================================================================
# MAIN PLOTTING
# ==========================================================================

# Detect iterations
iter_dirs = filter(readdir(output_dir)) do name
    startswith(name, "iteration_") && isdir(joinpath(output_dir, name))
end
iterations = sort([parse(Int, split(d, "_")[2]) for d in iter_dirs])

# Find last complete iteration (has G_ensemble.jld2)
last_complete = 0
for iter in reverse(iterations)
    iter_path = joinpath(output_dir, @sprintf("iteration_%03d", iter))
    if isfile(joinpath(iter_path, "G_ensemble.jld2"))
        last_complete = iter
        break
    end
end

if last_complete < 1
    @warn "Need at least iteration 1 to be complete for comparison. Using iteration 1 if available."
    last_complete = min(1, maximum(iterations))
end

initial_iter = 0
calibrated_iter = last_complete

@info "Comparing iteration $initial_iter (initial) vs iteration $calibrated_iter (calibrated)"

# Create figure
fig = Figure(size = (1500, 500 * length(short_names)))

for (row, short_name) in enumerate(short_names)
    @info "Processing $short_name..."

    # Load ERA5 observation
    era5_var = load_era5_var(short_name)

    # Load ensemble means
    initial_mean = load_ensemble_mean(output_dir, initial_iter, short_name)
    calibrated_mean = load_ensemble_mean(output_dir, calibrated_iter, short_name)

    # Resample to same grid (use ERA5 grid)
    era5_lons = ClimaAnalysis.longitudes(era5_var)
    era5_lats = ClimaAnalysis.latitudes(era5_var)

    initial_resampled = ClimaAnalysis.resampled_as(initial_mean; lon=era5_lons, lat=era5_lats)
    calibrated_resampled = ClimaAnalysis.resampled_as(calibrated_mean; lon=era5_lons, lat=era5_lats)

    # Set units to match ERA5 exactly (required for plot_bias_on_globe!)
    # Use the actual units from the ERA5 var to avoid string mismatch issues
    era5_units = ClimaAnalysis.units(era5_var)
    initial_resampled = ClimaAnalysis.set_units(initial_resampled, era5_units)
    calibrated_resampled = ClimaAnalysis.set_units(calibrated_resampled, era5_units)

    # Compute statistics
    initial_rmse = compute_global_rmse(initial_resampled, era5_var)
    initial_bias = compute_global_bias(initial_resampled, era5_var)
    calibrated_rmse = compute_global_rmse(calibrated_resampled, era5_var)
    calibrated_bias = compute_global_bias(calibrated_resampled, era5_var)

    @info "  Initial:    RMSE=$(round(initial_rmse, sigdigits=3)), Bias=$(round(initial_bias, sigdigits=3))"
    @info "  Calibrated: RMSE=$(round(calibrated_rmse, sigdigits=3)), Bias=$(round(calibrated_bias, sigdigits=3))"

    # Compute bias maps
    initial_bias_map = initial_resampled - era5_var
    calibrated_bias_map = calibrated_resampled - era5_var

    # Compute improvement map: |calibrated - obs| - |initial - obs|
    # Negative = improvement, Positive = worsened
    initial_abs_err = ClimaAnalysis.replace(abs, initial_bias_map)
    calibrated_abs_err = ClimaAnalysis.replace(abs, calibrated_bias_map)
    improvement_map = calibrated_abs_err - initial_abs_err

    # Get colorbar ranges
    bias_range = get(bias_ranges, short_name, extrema(initial_bias_map.data))
    improvement_range = get(improvement_ranges, short_name, extrema(improvement_map.data))

    # Panel 1: Initial bias - use plot_bias_on_globe! which handles bias calculation
    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        fig[row, 1],
        initial_resampled,
        era5_var;
        cmap_extrema = bias_range,
    )
    # Add custom title
    Label(fig[row, 1, Top()],
        "Initial (iter $initial_iter): $short_name\nRMSE=$(round(initial_rmse, sigdigits=3)), Bias=$(round(initial_bias, sigdigits=3))",
        fontsize = 12)

    # Panel 2: Calibrated bias
    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        fig[row, 2],
        calibrated_resampled,
        era5_var;
        cmap_extrema = bias_range,
    )
    Label(fig[row, 2, Top()],
        "Calibrated (iter $calibrated_iter): $short_name\nRMSE=$(round(calibrated_rmse, sigdigits=3)), Bias=$(round(calibrated_bias, sigdigits=3))",
        fontsize = 12)

    # Panel 3: Improvement - use heatmap2D_on_globe! with more_kwargs for custom colormap
    ClimaAnalysis.Visualize.heatmap2D_on_globe!(
        fig[row, 3],
        improvement_map;
        more_kwargs = Dict(:colormap => :RdBu, :colorrange => improvement_range),
    )
    Label(fig[row, 3, Top()],
        "Improvement: $short_name\n(blue=better, red=worse)",
        fontsize = 12)
end

# Save figure
output_path = joinpath(plots_dir, "rmse_comparison_$(exp_suffix).png")
save(output_path, fig)
@info "Saved: $output_path"

@info """
=====================================
RMSE Comparison Complete!
=====================================
Compared iteration $initial_iter (initial) vs iteration $calibrated_iter (calibrated)

Figure saved to: $output_path
"""
