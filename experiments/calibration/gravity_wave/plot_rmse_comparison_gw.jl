#=
plot_rmse_comparison_gw.jl

Compare initial vs calibrated ensemble bias against ERA5 observations
for 3D pressure-level variables (ta, ua, va) from gravity wave calibration.

Adapted from the plot_bias function in observation_map.jl.

Usage:
    julia --project=experiments/ClimaEarth plot_rmse_comparison_gw.jl [output_dir]

If output_dir is not provided, defaults to "output/gravity_wave"
=#

using ClimaAnalysis
using CairoMakie
using GeoMakie
using JLD2
using Statistics
using Printf
using DataStructures: OrderedDict
import Dates
import ClimaCoupler

# Override JLD2's default_iotype to avoid Lustre issues
JLD2.default_iotype() = IOStream

# ==========================================================================
# CONFIGURATION
# ==========================================================================
output_dir = length(ARGS) >= 1 ? ARGS[1] : "output/gravity_wave"
output_dir = abspath(output_dir)

# Create plots subdirectory
plots_dir = joinpath(output_dir, "plots")
mkpath(plots_dir)

exp_suffix = basename(output_dir)

@info "RMSE Comparison Plot (Gravity Wave Calibration)"
@info "Output directory: $output_dir"

# 3D pressure-level variables
short_names = ["ta", "ua", "va"]

# Pressure level to slice for 2D visualization (hPa)
# Set to nothing to auto-select the level with best improvement
PLOT_PRESSURE_LEVEL = nothing  # Will be determined per-variable

# ERA5 pressure levels used in calibration (must match observation_map.jl)
const ERA5_PRESSURE_LEVELS = [850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 250.0, 200.0]

# ERA5 uses 1979-01-01 as the time epoch
const ERA5_EPOCH = Dates.DateTime(1979, 1, 1)

# Units for each variable
var_units = Dict(
    "ta" => "K",
    "ua" => "m s^-1",
    "va" => "m s^-1",
)

# Colorbar ranges for bias plots
bias_ranges = Dict(
    "ta" => (-10, 10),
    "ua" => (-10, 10),
    "va" => (-5, 5),
)

# Colorbar ranges for improvement plots
improvement_ranges = Dict(
    "ta" => (-5, 5),
    "ua" => (-5, 5),
    "va" => (-2.5, 2.5),
)

# ==========================================================================
# DATA LOADING (adapted from observation_map.jl)
# ==========================================================================

"""
    load_era5_var(short_name)

Load preprocessed ERA5 variable from gravity_wave folder.
"""
function load_era5_var(short_name)
    preprocessed_path = joinpath(
        pkgdir(ClimaCoupler),
        "experiments/calibration/gravity_wave/preprocessed_vars.jld2"
    )

    if !isfile(preprocessed_path)
        error("Preprocessed vars not found at $preprocessed_path")
    end

    preprocessed_vars = JLD2.load_object(preprocessed_path)

    for var in preprocessed_vars
        if ClimaAnalysis.short_name(var) == short_name
            @info "  Loaded ERA5 $short_name"
            return var
        end
    end

    error("No ERA5 var found for $short_name")
end

"""
    get_var(short_name, simdir)

Load simulation variable and convert z-levels to pressure coordinates.
Adapted from observation_map.jl.
"""
function get_var(short_name, simdir)
    var = ClimaAnalysis.get(simdir; short_name, period="1M")

    # Convert 3D variables from z-levels to pressure coordinates
    if ClimaAnalysis.has_altitude(var)
        pfull = ClimaAnalysis.get(simdir; short_name="pfull", period="1M")
        # Window to avoid pressure inversions at low elevation
        pfull_windowed = ClimaAnalysis.window(pfull, "z", left=80)
        var_windowed = ClimaAnalysis.window(var, "z", left=80)
        var = ClimaAnalysis.Atmos.to_pressure_coordinates(var_windowed, pfull_windowed)
        var = ClimaAnalysis.Var.convert_dim_units(
            var,
            "pfull",
            "hPa";
            conversion_function = x -> 0.01 * x,
        )
    end

    var.attributes["short_name"] = short_name
    return var
end

"""
    preprocess_sim_var(var, short_name)

Preprocess simulation variable: set units, average time, resample to ERA5 pressure levels,
and rename dimensions to match ERA5 format.
Adapted from observation_map.jl plot_bias function.
"""
function preprocess_sim_var(var, short_name)
    # Set units
    var = ClimaAnalysis.set_units(var, var_units[short_name])

    # Handle time dimension
    if ClimaAnalysis.has_time(var)
        time_vals = ClimaAnalysis.times(var)
        if length(time_vals) == 1
            time_val = time_vals[1]
            var_3d = ClimaAnalysis.slice(var; time = time_val)
        else
            # Multiple time steps - average
            var_3d = ClimaAnalysis.average_time(var)
        end
    else
        var_3d = var
    end

    # Resample to ERA5 pressure levels
    if ClimaAnalysis.has_pressure(var_3d)
        lon = ClimaAnalysis.longitudes(var_3d)
        lat = ClimaAnalysis.latitudes(var_3d)
        var_resampled = ClimaAnalysis.resampled_as(var_3d; lon, lat, pressure_level = ERA5_PRESSURE_LEVELS)

        # Rename dimensions to match ERA5 observation format
        renamed_dims = OrderedDict{String, AbstractArray}()
        for (k, v) in var_resampled.dims
            new_key = if k == "pfull"
                "pressure_level"
            elseif k == "lon"
                "longitude"
            elseif k == "lat"
                "latitude"
            else
                k
            end
            renamed_dims[new_key] = v
        end

        renamed_dim_attribs = OrderedDict{String, AbstractDict}()
        for (k, v) in var_resampled.dim_attributes
            new_key = if k == "pfull"
                "pressure_level"
            elseif k == "lon"
                "longitude"
            elseif k == "lat"
                "latitude"
            else
                k
            end
            renamed_dim_attribs[new_key] = v
        end

        var_out = ClimaAnalysis.OutputVar(renamed_dims, var_resampled.data)
        merge!(var_out.attributes, var_resampled.attributes)
        merge!(var_out.dim_attributes, renamed_dim_attribs)
        var_out.attributes["short_name"] = short_name
        return var_out
    else
        return var_3d
    end
end

"""
    load_ensemble_mean(output_dir, iteration, short_name)

Load simulation outputs for all members at given iteration and compute ensemble mean.
"""
function load_ensemble_mean(output_dir, iteration, short_name)
    iter_str = @sprintf("iteration_%03d", iteration)
    iter_path = joinpath(output_dir, iter_str)

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
        simdir_path = joinpath(member_path, "amip_land", "output_active")

        if !isdir(simdir_path)
            @warn "No output_active directory at $simdir_path, skipping"
            continue
        end

        try
            simdir = ClimaAnalysis.SimDir(simdir_path)
            var = get_var(short_name, simdir)
            var = preprocess_sim_var(var, short_name)
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
    ensemble_mean.attributes["short_name"] = short_name

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

    # Load ERA5 observation (keep 3D for level selection)
    era5_var_3d = load_era5_var(short_name)

    # Average over time if present
    if ClimaAnalysis.has_time(era5_var_3d)
        n_times = length(ClimaAnalysis.times(era5_var_3d))
        if n_times > 1
            era5_var_3d = ClimaAnalysis.average_time(era5_var_3d)
        else
            era5_var_3d = ClimaAnalysis.slice(era5_var_3d; time = first(ClimaAnalysis.times(era5_var_3d)))
        end
    end

    # Load ensemble means (keep 3D for level selection)
    initial_mean_3d = load_ensemble_mean(output_dir, initial_iter, short_name)
    calibrated_mean_3d = load_ensemble_mean(output_dir, calibrated_iter, short_name)

    # Find the best pressure level (most improvement or least worsening)
    @info "  Scanning pressure levels for $short_name..."
    best_level = ERA5_PRESSURE_LEVELS[1]
    best_improvement = Inf  # Lower is better (negative = improvement)

    for plev in ERA5_PRESSURE_LEVELS
        # Slice at this level
        era5_at_lev = ClimaAnalysis.slice(era5_var_3d; pressure_level = plev)
        init_at_lev = ClimaAnalysis.slice(initial_mean_3d; pressure_level = plev)
        calib_at_lev = ClimaAnalysis.slice(calibrated_mean_3d; pressure_level = plev)

        # Resample to ERA5 grid
        era5_lons = ClimaAnalysis.longitudes(era5_at_lev)
        era5_lats = ClimaAnalysis.latitudes(era5_at_lev)
        init_resampled = ClimaAnalysis.resampled_as(init_at_lev; longitude=era5_lons, latitude=era5_lats)
        calib_resampled = ClimaAnalysis.resampled_as(calib_at_lev; longitude=era5_lons, latitude=era5_lats)

        # Set units
        era5_units = ClimaAnalysis.units(era5_at_lev)
        init_resampled = ClimaAnalysis.set_units(init_resampled, era5_units)
        calib_resampled = ClimaAnalysis.set_units(calib_resampled, era5_units)

        # Compute RMSE
        init_rmse = compute_global_rmse(init_resampled, era5_at_lev)
        calib_rmse = compute_global_rmse(calib_resampled, era5_at_lev)
        delta_rmse = calib_rmse - init_rmse  # Negative = improvement

        improvement_pct = 100 * delta_rmse / init_rmse
        status = delta_rmse < 0 ? "✓" : "✗"
        @info @sprintf("    %4.0f hPa: RMSE %.3f → %.3f (Δ = %+.3f, %+.1f%%) %s",
                       plev, init_rmse, calib_rmse, delta_rmse, improvement_pct, status)

        if delta_rmse < best_improvement
            best_improvement = delta_rmse
            best_level = plev
        end
    end

    # Use specified level or best level
    plot_level = isnothing(PLOT_PRESSURE_LEVEL) ? best_level : PLOT_PRESSURE_LEVEL
    @info "  Selected pressure level: $(Int(plot_level)) hPa (Δ RMSE = $(round(best_improvement, sigdigits=3)))"

    # Slice at plotting pressure level
    era5_var = ClimaAnalysis.slice(era5_var_3d; pressure_level = plot_level)
    initial_mean = ClimaAnalysis.slice(initial_mean_3d; pressure_level = plot_level)
    calibrated_mean = ClimaAnalysis.slice(calibrated_mean_3d; pressure_level = plot_level)

    # Resample to same grid (use ERA5 grid)
    era5_lons = ClimaAnalysis.longitudes(era5_var)
    era5_lats = ClimaAnalysis.latitudes(era5_var)

    # Note: After preprocess_sim_var, dims are named longitude/latitude
    initial_resampled = ClimaAnalysis.resampled_as(initial_mean; longitude=era5_lons, latitude=era5_lats)
    calibrated_resampled = ClimaAnalysis.resampled_as(calibrated_mean; longitude=era5_lons, latitude=era5_lats)

    # Set units to match ERA5 exactly
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

    # Panel 1: Initial bias
    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        fig[row, 1],
        initial_resampled,
        era5_var;
        cmap_extrema = bias_range,
    )
    Label(fig[row, 1, Top()],
        "Initial (iter $initial_iter): $short_name @ $(Int(plot_level)) hPa\nRMSE=$(round(initial_rmse, sigdigits=3)), Bias=$(round(initial_bias, sigdigits=3))",
        fontsize = 12)

    # Panel 2: Calibrated bias
    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        fig[row, 2],
        calibrated_resampled,
        era5_var;
        cmap_extrema = bias_range,
    )
    Label(fig[row, 2, Top()],
        "Calibrated (iter $calibrated_iter): $short_name @ $(Int(plot_level)) hPa\nRMSE=$(round(calibrated_rmse, sigdigits=3)), Bias=$(round(calibrated_bias, sigdigits=3))",
        fontsize = 12)

    # Panel 3: Improvement
    ClimaAnalysis.Visualize.heatmap2D_on_globe!(
        fig[row, 3],
        improvement_map;
        more_kwargs = Dict(:colormap => :RdBu, :colorrange => improvement_range),
    )
    Label(fig[row, 3, Top()],
        "Improvement: $short_name @ $(Int(plot_level)) hPa\n(blue=better, red=worse)",
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
Pressure levels: auto-selected per variable (best improvement)

Figure saved to: $output_path
"""
