#=
plot_bias_comparison.jl

Compare bias plots between initial and calibrated iterations side-by-side,
using the existing infrastructure from observation_map.jl.

Layout: Each row is a variable, columns are [Initial, Calibrated, Improvement]

Usage:
    julia --project=experiments/ClimaEarth plot_bias_comparison.jl [output_dir] [calibrated_iter] [initial_iter]

    output_dir:      Path to calibration output (default: output/gravity_wave)
    calibrated_iter: Iteration to use as "calibrated" (default: auto-detect last complete)
    initial_iter:    Iteration to use as "initial" (default: 0)
=#

using Dates
using DataStructures: OrderedDict
using Statistics: quantile
import ClimaAnalysis
import ClimaAnalysis: short_name, slice, window, latitudes, longitudes
import ClimaAnalysis: set_units, resampled_as, weighted_average_lonlat
import ClimaCalibrate
import ClimaCoupler
import EnsembleKalmanProcesses as EKP
import GeoMakie
import JLD2
using Printf
using Statistics

# Load the calibration config and observation_map infrastructure
include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "api.jl"))
include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "gravity_wave", "run_calibration.jl"))
include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "gravity_wave", "observation_map.jl"))

"""
    format_units(units_str)

Convert unit strings to Unicode for better display (e.g., "m s^-1" → "m s⁻¹").
"""
function format_units(units_str)
    # Unicode superscript mappings
    superscripts = Dict(
        '^' => "",  # Remove the caret
        '-' => '⁻',
        '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
        '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
    )

    result = ""
    in_superscript = false
    for c in units_str
        if c == '^'
            in_superscript = true
        elseif in_superscript && haskey(superscripts, c)
            result *= superscripts[c]
        elseif c == ' ' || c == '*'
            in_superscript = false
            result *= c
        else
            in_superscript = false
            result *= c
        end
    end
    return result
end

# Override JLD2's default_iotype to avoid Lustre issues
JLD2.default_iotype() = IOStream

# ==========================================================================
# CONFIGURATION
# ==========================================================================
output_dir = length(ARGS) >= 1 ? ARGS[1] : "output/gravity_wave"
output_dir = abspath(output_dir)

# Optional: specify calibrated and initial iterations explicitly
calibrated_iter_arg = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nothing
initial_iter_arg = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : nothing

plots_dir = joinpath(output_dir, "plots")
mkpath(plots_dir)

exp_suffix = basename(output_dir)

# Pressure level for 2D slicing (hPa)
# Set to nothing to auto-select the level with best improvement per variable
PLOT_PRESSURE_LEVEL = nothing

# Include va (meridional wind) in plots? Set to false for GW-focused analysis.
INCLUDE_VA = false

@info "Bias Comparison Plot"
@info "Output directory: $output_dir"
@info "Include va: $INCLUDE_VA"

# ==========================================================================
# AUTO-DETECT ITERATIONS
# ==========================================================================
function detect_iterations(output_dir; skip_latest=true)
    iter_dirs = filter(readdir(output_dir)) do name
        startswith(name, "iteration_") && isdir(joinpath(output_dir, name))
    end
    iterations = sort([parse(Int, split(d, "_")[2]) for d in iter_dirs])

    if isempty(iterations)
        error("No iteration directories found in $output_dir")
    end

    # Find iterations that have eki_file.jld2 (indicates EKP update completed)
    complete_iters = filter(iterations) do iter
        iter_path = joinpath(output_dir, @sprintf("iteration_%03d", iter))
        isfile(joinpath(iter_path, "eki_file.jld2"))
    end

    if isempty(complete_iters)
        @warn "No iterations with eki_file.jld2 found, using iteration 0"
        return iterations, 0
    end

    # If skip_latest is true and we have more than one complete iteration,
    # skip the latest one (it might still be running)
    if skip_latest && length(complete_iters) > 1
        last_complete = complete_iters[end-1]
        @info "Skipping latest iteration $(complete_iters[end]) (may still be running)"
    else
        last_complete = complete_iters[end]
    end

    return iterations, last_complete
end

# ==========================================================================
# HELPER: Get sim and obs vars for an iteration (adapted from plot_bias)
# ==========================================================================
function get_vars_for_iteration(output_dir, iteration)
    iter_str = @sprintf("iteration_%03d", iteration)
    ekp_path = joinpath(output_dir, iter_str, "eki_file.jld2")
    ekp = JLD2.load(ekp_path)["single_stored_object"]

    # Get observations from EKP
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges[iteration + 1]
    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs = ClimaCalibrate.ObservationRecipe.get_observations_for_nth_iteration(
        obs_series,
        iteration + 1,
    )

    # Reconstruct OutputVars from observations (ERA5 data)
    era5_vars = []
    for obs in minibatch_obs
        obs_vars = ClimaCalibrate.ObservationRecipe.reconstruct_vars(obs)
        append!(era5_vars, obs_vars)
    end

    # Find first member's simdir
    iter_path = joinpath(output_dir, iter_str)
    member_dirs = filter(readdir(iter_path)) do name
        startswith(name, "member_") && isdir(joinpath(iter_path, name))
    end
    member_path = joinpath(iter_path, first(sort(member_dirs)), "amip_land/output_active")
    simdir = ClimaAnalysis.SimDir(member_path)

    # Get simulation variables (same preprocessing as plot_bias)
    # Filter short_names based on INCLUDE_VA setting
    plot_short_names = INCLUDE_VA ? CALIBRATE_CONFIG.short_names : filter(s -> s != "va", CALIBRATE_CONFIG.short_names)
    sim_vars = []
    for sn in plot_short_names
        var = get_var(sn, simdir)
        period = largest_period(sample_date_ranges)
        var = ClimaAnalysis.Var._shift_by(var, date -> date - period)
        var = set_units(var, var_units[sn])
        var = window(var, "time"; left = sample_date_ranges[1], right = sample_date_ranges[2])

        # Resample 3D pressure-level variables to ERA5 standard levels
        if ClimaAnalysis.has_pressure(var)
            if ClimaAnalysis.has_time(var)
                time_vals = ClimaAnalysis.times(var)
                if length(time_vals) == 1
                    time_val = time_vals[1]
                    var_3d = slice(var; time = time_val)
                    lon = longitudes(var_3d)
                    lat = latitudes(var_3d)
                    var_resampled = resampled_as(var_3d; lon, lat, pressure_level = ERA5_PRESSURE_LEVELS)

                    # Rename dimensions to match ERA5
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

                    # Convert time
                    sim_start_date_str = get(var.attributes, "start_date", nothing)
                    if !isnothing(sim_start_date_str)
                        sim_start_date = Dates.DateTime(sim_start_date_str)
                        absolute_datetime = sim_start_date + Dates.Second(round(Int, time_val))
                        era5_time_val = Float64(Dates.value(absolute_datetime - ERA5_EPOCH) / 1000)
                    else
                        era5_time_val = time_val
                    end

                    time_dim = OrderedDict("time" => [era5_time_val])
                    new_dims = merge(renamed_dims, time_dim)
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
                    time_dim_attribs = OrderedDict("time" => Dict{String, Any}("units" => "s"))
                    new_dim_attribs = merge(renamed_dim_attribs, time_dim_attribs)
                    new_data = reshape(var_resampled.data, size(var_resampled.data)..., 1)
                    var = ClimaAnalysis.OutputVar(new_dims, new_data)
                    merge!(var.attributes, var_resampled.attributes)
                    merge!(var.dim_attributes, new_dim_attribs)
                    var.attributes["start_date"] = string(ERA5_EPOCH)
                end
            end
        end
        push!(sim_vars, var)
    end

    # Match sim with era5 by short_name
    # Return 3D vars (without pressure slicing) for level selection later
    var_pairs = []
    sample_date = first(unique(sample_date_ranges))
    for sim_var in sim_vars
        sn = short_name(sim_var)
        era5_idx = findfirst(v -> short_name(v) == sn, era5_vars)
        if !isnothing(era5_idx)
            sim_t = slice(sim_var; time = sample_date)
            era5_t = slice(era5_vars[era5_idx]; time = sample_date)
            # Don't slice at pressure level here - return 3D for level selection
            push!(var_pairs, (sn, sim_t, era5_t))
        end
    end

    return var_pairs
end

"""
    compute_rmse(sim_var, obs_var)

Compute global RMSE (weighted by latitude).
"""
function compute_rmse(sim_var, obs_var)
    # Resample sim to obs grid
    obs_lons = ClimaAnalysis.longitudes(obs_var)
    obs_lats = ClimaAnalysis.latitudes(obs_var)
    sim_r = ClimaAnalysis.resampled_as(sim_var; longitude = obs_lons, latitude = obs_lats)

    diff = sim_r - obs_var
    sq_diff = ClimaAnalysis.replace(x -> x^2, diff)
    mean_sq = ClimaAnalysis.average_lonlat(sq_diff; weighted = true)
    return sqrt(mean_sq.data[1])
end

"""
    find_best_pressure_level(init_sim, init_obs, calib_sim, calib_obs)

Scan pressure levels and find the one with best improvement (most negative Δ RMSE).
Returns (best_level, best_improvement).
"""
function find_best_pressure_level(init_sim, init_obs, calib_sim, calib_obs)
    best_level = ERA5_PRESSURE_LEVELS[1]
    best_improvement = Inf

    for plev in ERA5_PRESSURE_LEVELS
        init_sim_lev = slice(init_sim; pressure_level = plev)
        init_obs_lev = slice(init_obs; pressure_level = plev)
        calib_sim_lev = slice(calib_sim; pressure_level = plev)
        calib_obs_lev = slice(calib_obs; pressure_level = plev)

        init_rmse = compute_rmse(init_sim_lev, init_obs_lev)
        calib_rmse = compute_rmse(calib_sim_lev, calib_obs_lev)
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

    return best_level, best_improvement
end

# ==========================================================================
# MAIN
# ==========================================================================
# Use explicit iteration if provided, otherwise auto-detect
skip_latest = isnothing(calibrated_iter_arg)
iterations, last_complete = detect_iterations(output_dir; skip_latest)
@info "Found iterations: $iterations, last complete (auto): $last_complete"

initial_iter = isnothing(initial_iter_arg) ? 0 : initial_iter_arg
calibrated_iter = isnothing(calibrated_iter_arg) ? last_complete : calibrated_iter_arg

@info "Comparing iteration $initial_iter (initial) vs iteration $calibrated_iter (calibrated)"

# Get vars for each iteration
@info "Loading initial iteration data..."
initial_pairs = get_vars_for_iteration(output_dir, initial_iter)

@info "Loading calibrated iteration data..."
calibrated_pairs = get_vars_for_iteration(output_dir, calibrated_iter)

# Create comparison figure (smaller = text appears larger)
n_vars = length(initial_pairs)
# figure_padding = (left, right, bottom, top)
fig = GeoMakie.Figure(size = (1200, 80 + 230 * n_vars); figure_padding = (5, 5, 40, 50))

# Column titles (use row 1 position with offset, will adjust with rowsize)
GeoMakie.Label(fig[0, 1], "Initial (iter $initial_iter)", fontsize = 24, font = :bold)
GeoMakie.Label(fig[0, 2], "Calibrated (iter $calibrated_iter)", fontsize = 24, font = :bold)
GeoMakie.Label(fig[0, 3], "Improvement\n(blue=better, red=worse)", fontsize = 24, font = :bold)

for (row, (sn, init_sim_3d, init_era5_3d)) in enumerate(initial_pairs)
    @info "Processing $sn..."

    # Find matching calibrated pair
    calib_idx = findfirst(p -> p[1] == sn, calibrated_pairs)
    if isnothing(calib_idx)
        @warn "No calibrated data for $sn"
        continue
    end
    _, calib_sim_3d, calib_era5_3d = calibrated_pairs[calib_idx]

    # Determine pressure level to use
    if ClimaAnalysis.has_pressure(init_sim_3d)
        if isnothing(PLOT_PRESSURE_LEVEL)
            # Scan pressure levels and find the best one
            @info "  Scanning pressure levels for $sn..."
            plot_level, delta_rmse = find_best_pressure_level(
                init_sim_3d, init_era5_3d, calib_sim_3d, calib_era5_3d
            )
            @info "  Selected: $(Int(plot_level)) hPa (Δ RMSE = $(round(delta_rmse, sigdigits=3)))"
        else
            plot_level = PLOT_PRESSURE_LEVEL
        end

        # Slice at selected pressure level
        init_sim = slice(init_sim_3d; pressure_level = plot_level)
        init_era5 = slice(init_era5_3d; pressure_level = plot_level)
        calib_sim = slice(calib_sim_3d; pressure_level = plot_level)
        calib_era5 = slice(calib_era5_3d; pressure_level = plot_level)
        level_str = " @ $(Int(plot_level)) hPa"
    else
        # 2D variable, no pressure slicing needed
        init_sim = init_sim_3d
        init_era5 = init_era5_3d
        calib_sim = calib_sim_3d
        calib_era5 = calib_era5_3d
        plot_level = nothing
        level_str = ""
    end

    # Compute biases
    init_bias = ClimaAnalysis.global_bias(init_sim, init_era5)
    calib_bias = ClimaAnalysis.global_bias(calib_sim, calib_era5)
    pct_change = 100 * (calib_bias - init_bias) / abs(init_bias)

    # Resample all vars to ERA5 grid for arithmetic compatibility
    era5_lons = ClimaAnalysis.longitudes(init_era5)
    era5_lats = ClimaAnalysis.latitudes(init_era5)

    init_sim_r = ClimaAnalysis.resampled_as(init_sim; longitude = era5_lons, latitude = era5_lats)
    calib_sim_r = ClimaAnalysis.resampled_as(calib_sim; longitude = era5_lons, latitude = era5_lats)

    # Compute bias maps (sim - obs)
    init_bias_map = init_sim_r - init_era5
    calib_bias_map = calib_sim_r - calib_era5

    # Get units
    units = ClimaAnalysis.units(init_sim)

    # Compute improvement: |calib - obs| - |init - obs|  (negative = improvement)
    init_abs_err = ClimaAnalysis.replace(abs, init_bias_map)
    calib_abs_err = ClimaAnalysis.replace(abs, calib_bias_map)
    improvement = calib_abs_err - init_abs_err

    # Compute colorbar range from initial bias
    # Row 1: Use 5th-95th percentiles to avoid Antarctic extremes dominating
    # Other rows: Use full range (min/max)
    data_flat = vec(init_bias_map.data)
    if row == 1
        q_low = quantile(data_flat, 0.05)
        q_high = quantile(data_flat, 0.95)
        q_max = max(abs(q_low), abs(q_high))
    else
        q_max = maximum(abs.(data_flat))
    end
    cmap_extrema = (-q_max, q_max)

    # Use GeoMakie directly for more control over the plot
    # Common axis options: hide lon/lat labels
    axis_opts = (; dest = "+proj=wintri", xticklabelsvisible = false, yticklabelsvisible = false)

    # Colormap: Reverse(:RdBu) so high values (warm/fast) are red
    cmap = GeoMakie.Reverse(:RdBu)

    # Panel 1: Initial bias (rasterize=2 for PDF compression - converts to bitmap)
    ax1 = GeoMakie.GeoAxis(fig[row, 1]; axis_opts...)
    lons = ClimaAnalysis.longitudes(init_bias_map)
    lats = ClimaAnalysis.latitudes(init_bias_map)
    data1 = init_bias_map.data
    GeoMakie.surface!(ax1, lons, lats, data1; colormap = cmap, colorrange = cmap_extrema, shading = GeoMakie.NoShading, rasterize = 2)
    GeoMakie.lines!(ax1, GeoMakie.coastlines(); color = :black, linewidth = 0.5)

    # Panel 2: Calibrated bias
    ax2 = GeoMakie.GeoAxis(fig[row, 2]; axis_opts...)
    data2 = calib_bias_map.data
    GeoMakie.surface!(ax2, lons, lats, data2; colormap = cmap, colorrange = cmap_extrema, shading = GeoMakie.NoShading, rasterize = 2)
    GeoMakie.lines!(ax2, GeoMakie.coastlines(); color = :black, linewidth = 0.5)

    # Panel 3: Improvement map (same colorrange as bias panels)
    ax3 = GeoMakie.GeoAxis(fig[row, 3]; axis_opts...)
    data3 = improvement.data
    GeoMakie.surface!(ax3, lons, lats, data3; colormap = cmap, colorrange = cmap_extrema, shading = GeoMakie.NoShading, rasterize = 2)
    GeoMakie.lines!(ax3, GeoMakie.coastlines(); color = :black, linewidth = 0.5)

    # Format units for display
    units_display = format_units(units)

    # Row label on left side (units shown on colorbar, not here)
    GeoMakie.Label(fig[row, 0],
        "$sn$level_str",
        fontsize = 24,
        rotation = π/2)

    # Single shared colorbar for all panels in the row
    GeoMakie.Colorbar(fig[row, 4],
        colormap = cmap,
        colorrange = cmap_extrema,
        label = "[$units_display]",
        width = 15,
        ticklabelsize = 24,
        labelsize = 24,
        labelpadding = -5)

    # Add bias values as text annotations
    GeoMakie.Label(fig[row, 1, GeoMakie.Bottom()],
        "Bias = $(round(init_bias, sigdigits=3))",
        fontsize = 24)
    GeoMakie.Label(fig[row, 2, GeoMakie.Bottom()],
        "Bias = $(round(calib_bias, sigdigits=3))",
        fontsize = 24)
    GeoMakie.Label(fig[row, 3, GeoMakie.Bottom()],
        "$(@sprintf("%+.1f%%", pct_change)) bias change",
        fontsize = 24)

    @info @sprintf("  %s%s: Bias %.3f → %.3f (%+.1f%%)", sn, level_str, init_bias, calib_bias, pct_change)
end

# Set layout sizes AFTER content is placed
GeoMakie.colsize!(fig.layout, 0, GeoMakie.Fixed(55))   # Row labels
GeoMakie.colsize!(fig.layout, 1, GeoMakie.Fixed(370))  # Panel 1
GeoMakie.colsize!(fig.layout, 2, GeoMakie.Fixed(370))  # Panel 2
GeoMakie.colsize!(fig.layout, 3, GeoMakie.Fixed(370))  # Panel 3
GeoMakie.colsize!(fig.layout, 4, GeoMakie.Fixed(65))   # Colorbar

GeoMakie.rowsize!(fig.layout, 0, GeoMakie.Fixed(45))   # Title row
for r in 1:n_vars
    GeoMakie.rowsize!(fig.layout, r, GeoMakie.Fixed(230))  # Data rows
end

GeoMakie.colgap!(fig.layout, -5)  # Negative gap to overlap slightly
GeoMakie.rowgap!(fig.layout, -5)  # Negative gap to reduce space

# Resize figure to fit all content (prevents clipping)
GeoMakie.resize_to_layout!(fig)

# Save as PNG and PDF (rasterize=2 on surface! makes PDF small)
output_path_png = joinpath(plots_dir, "bias_comparison_$(exp_suffix).png")
output_path_pdf = joinpath(plots_dir, "bias_comparison_$(exp_suffix).pdf")
GeoMakie.save(output_path_png, fig)
GeoMakie.save(output_path_pdf, fig)
@info "Saved: $output_path_png"
@info "Saved: $output_path_pdf"

level_info = isnothing(PLOT_PRESSURE_LEVEL) ? "auto-selected per variable (best improvement)" : "$(Int(PLOT_PRESSURE_LEVEL)) hPa"
@info """
=====================================
Bias Comparison Complete!
=====================================
Compared iteration $initial_iter (initial) vs iteration $calibrated_iter (calibrated)
Pressure levels: $level_info

Figures saved to:
  $output_path_png
  $output_path_pdf
"""
