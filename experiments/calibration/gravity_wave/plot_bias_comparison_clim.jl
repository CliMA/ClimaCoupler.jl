#=
plot_bias_comparison_clim.jl

Show bias maps for 2010, 2011, 2012 (ensemble means) vs ERA5 February climatology.

Each row is a variable (ta, ua), each column is a year (2010, 2011, 2012).
All compared against the same multi-year February climatology reference.

Usage:
    julia --project=experiments/ClimaEarth plot_bias_comparison_clim.jl [output_dir]

    output_dir: Path to calibration output (default: output/gravity_wave)
=#

using Dates
using DataStructures: OrderedDict
using LaTeXStrings
using Statistics: quantile
import ClimaAnalysis
import ClimaAnalysis: short_name, slice, window, latitudes, longitudes
import ClimaAnalysis: set_units, resampled_as
import ClimaCoupler
import GeoMakie
import JLD2
using Printf
using Statistics

# Load the calibration config and observation_map infrastructure
include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "api.jl"))
include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "gravity_wave", "run_calibration.jl"))
include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "gravity_wave", "observation_map.jl"))

# LaTeX-formatted units for colorbar labels
const UNITS_LATEX = Dict(
    "ta" => L"K",
    "ua" => L"m s$^{-1}$",
    "va" => L"m s$^{-1}$",
)

# Override JLD2's default_iotype to avoid Lustre issues
JLD2.default_iotype() = IOStream

# ==========================================================================
# CONFIGURATION
# ==========================================================================
output_dir = length(ARGS) >= 1 ? ARGS[1] : "output/gravity_wave"
output_dir = abspath(output_dir)

plots_dir = joinpath(output_dir, "plots")
mkpath(plots_dir)

exp_suffix = basename(output_dir)

# Years/iterations to plot
years = [2010, 2011, 2012]
iterations = [0, 1, 2]

# Pressure level for 2D slicing (hPa)
PLOT_PRESSURE_LEVEL = 500.0

# Variables to plot (exclude va for GW-focused analysis)
plot_short_names = ["ta", "ua"]

@info "Bias Comparison Plot (using February Climatology)"
@info "Output directory: $output_dir"
@info "Pressure level: $(Int(PLOT_PRESSURE_LEVEL)) hPa"

# ==========================================================================
# ERA5 CLIMATOLOGY LOADING
# ==========================================================================

const ERA5_EPOCH_CLIM = Dates.DateTime(1979, 1, 1)

"""
Load ERA5 observation variable and compute February climatology from all available years.
"""
function load_era5_climatology(short_name, target_month=2)
    preprocessed_path = joinpath(
        pkgdir(ClimaCoupler),
        "experiments/calibration/gravity_wave/preprocessed_vars.jld2"
    )

    if !isfile(preprocessed_path)
        error("Preprocessed vars not found at $preprocessed_path")
    end

    preprocessed_vars = JLD2.load_object(preprocessed_path)

    era5_var = nothing
    for var in preprocessed_vars
        if ClimaAnalysis.short_name(var) == short_name
            era5_var = var
            break
        end
    end

    if isnothing(era5_var)
        error("No ERA5 var found for $short_name")
    end

    if !ClimaAnalysis.has_time(era5_var)
        return era5_var
    end

    # Find all time points that correspond to the target month
    time_vals = ClimaAnalysis.times(era5_var)
    feb_indices = Int[]

    for (i, t) in enumerate(time_vals)
        dt = ERA5_EPOCH_CLIM + Dates.Second(round(Int, t))
        if Dates.month(dt) == target_month
            push!(feb_indices, i)
        end
    end

    if isempty(feb_indices)
        @warn "No $(Dates.monthname(target_month)) data found, using time average"
        return ClimaAnalysis.average_time(era5_var)
    end

    @info "  Computing $(Dates.monthname(target_month)) climatology from $(length(feb_indices)) years"

    # Average all February months to get climatology
    feb_vars = [ClimaAnalysis.slice(era5_var; time = time_vals[i]) for i in feb_indices]

    climatology = feb_vars[1]
    for i in 2:length(feb_vars)
        climatology = climatology + feb_vars[i]
    end
    climatology = climatology / length(feb_vars)
    climatology.attributes["short_name"] = short_name

    # Preserve units from original variable
    original_units = ClimaAnalysis.units(era5_var)
    if !isempty(original_units)
        climatology = ClimaAnalysis.set_units(climatology, original_units)
    end

    return climatology
end

# ==========================================================================
# ENSEMBLE MEAN LOADING (adapted from plot_rmse_barchart.jl)
# ==========================================================================

"""
Load simulation variable and convert z-levels to pressure coordinates.
"""
function get_var_local(short_name, simdir)
    var = ClimaAnalysis.get(simdir; short_name, period="1M")

    if ClimaAnalysis.has_altitude(var)
        pfull = ClimaAnalysis.get(simdir; short_name="pfull", period="1M")
        pfull_windowed = ClimaAnalysis.window(pfull, "z", left=80)
        var_windowed = ClimaAnalysis.window(var, "z", left=80)
        var = ClimaAnalysis.Atmos.to_pressure_coordinates(var_windowed, pfull_windowed)
        var = ClimaAnalysis.Var.convert_dim_units(
            var, "pfull", "hPa"; conversion_function = x -> 0.01 * x
        )
    end

    var.attributes["short_name"] = short_name
    return var
end

"""
Preprocess simulation variable for comparison with ERA5.
"""
function preprocess_sim_var(var, short_name)
    var = ClimaAnalysis.set_units(var, var_units[short_name])

    # Handle time dimension
    if ClimaAnalysis.has_time(var)
        time_vals = ClimaAnalysis.times(var)
        if length(time_vals) == 1
            var_3d = ClimaAnalysis.slice(var; time = time_vals[1])
        else
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

        # Rename dimensions to match ERA5
        renamed_dims = OrderedDict{String, AbstractArray}()
        for (k, v) in var_resampled.dims
            new_key = k == "pfull" ? "pressure_level" : (k == "lon" ? "longitude" : (k == "lat" ? "latitude" : k))
            renamed_dims[new_key] = v
        end

        renamed_dim_attribs = OrderedDict{String, AbstractDict}()
        for (k, v) in var_resampled.dim_attributes
            new_key = k == "pfull" ? "pressure_level" : (k == "lon" ? "longitude" : (k == "lat" ? "latitude" : k))
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
Load ensemble mean for a given iteration.
"""
function load_ensemble_mean_local(output_dir, iteration, short_name)
    iter_str = @sprintf("iteration_%03d", iteration)
    iter_path = joinpath(output_dir, iter_str)

    member_dirs = filter(readdir(iter_path)) do name
        startswith(name, "member_") && isdir(joinpath(iter_path, name))
    end

    if isempty(member_dirs)
        error("No member directories found in $iter_path")
    end

    member_vars = []
    for member_dir in sort(member_dirs)
        member_path = joinpath(iter_path, member_dir)
        simdir_path = joinpath(member_path, "amip_land", "output_active")

        if !isdir(simdir_path)
            continue
        end

        try
            simdir = ClimaAnalysis.SimDir(simdir_path)
            var = get_var_local(short_name, simdir)
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

    @info "  Loaded ensemble mean from $(length(member_vars)) members"
    return ensemble_mean
end

# ==========================================================================
# LOAD DATA
# ==========================================================================

@info "Loading ERA5 February climatology..."
era5_climatology = Dict{String, Any}()
for sn in plot_short_names
    era5_climatology[sn] = load_era5_climatology(sn, 2)
end

@info "Loading ensemble means for each iteration..."
sim_data = Dict{String, Vector{Any}}()  # sim_data[sn][iter_idx] = var
for sn in plot_short_names
    sim_data[sn] = []
    for (i, iteration) in enumerate(iterations)
        @info "  Loading $sn for iteration $iteration (year $(years[i]))..."
        var = load_ensemble_mean_local(output_dir, iteration, sn)
        push!(sim_data[sn], var)
    end
end

# ==========================================================================
# CREATE FIGURE: rows = variables, columns = years
# ==========================================================================

n_vars = length(plot_short_names)
n_years = length(years)

fig = GeoMakie.Figure(size = (1200, 80 + 230 * n_vars); figure_padding = (5, 5, 40, 50))

# Column titles
for (col, year) in enumerate(years)
    iter_label = col == 1 ? "(Initial)" : "(Calibrated)"
    GeoMakie.Label(fig[0, col], "$year\n$iter_label", fontsize = 24, font = :bold)
end

for (row, sn) in enumerate(plot_short_names)
    @info "Processing $sn..."

    era5_3d = era5_climatology[sn]

    # Slice at pressure level
    era5_ref = slice(era5_3d; pressure_level = PLOT_PRESSURE_LEVEL)
    era5_lons = ClimaAnalysis.longitudes(era5_ref)
    era5_lats = ClimaAnalysis.latitudes(era5_ref)

    # Compute bias maps for all years
    bias_maps = []
    biases = Float64[]

    for (col, iteration) in enumerate(iterations)
        sim_3d = sim_data[sn][col]
        sim_slice = slice(sim_3d; pressure_level = PLOT_PRESSURE_LEVEL)

        # Resample to ERA5 grid and match units
        sim_r = ClimaAnalysis.resampled_as(sim_slice; longitude = era5_lons, latitude = era5_lats)
        era5_units = ClimaAnalysis.units(era5_ref)
        sim_r = ClimaAnalysis.set_units(sim_r, era5_units)

        # Compute bias map
        bias_map = sim_r - era5_ref
        push!(bias_maps, bias_map)

        # Compute global bias
        bias_val = ClimaAnalysis.global_bias(sim_r, era5_ref)
        push!(biases, bias_val)
    end

    # Compute colorbar range from all bias maps (consistent across row)
    all_data = vcat([vec(bm.data) for bm in bias_maps]...)
    q_low = quantile(all_data, 0.05)
    q_high = quantile(all_data, 0.95)
    q_max = max(abs(q_low), abs(q_high))
    cmap_extrema = (-q_max, q_max)

    # Get LaTeX-formatted units for colorbar
    units_label = UNITS_LATEX[sn]

    # Axis options
    axis_opts = (; dest = "+proj=wintri", xticklabelsvisible = false, yticklabelsvisible = false)
    cmap = GeoMakie.Reverse(:RdBu)

    # Plot each year
    for (col, year) in enumerate(years)
        ax = GeoMakie.GeoAxis(fig[row, col]; axis_opts...)
        lons = ClimaAnalysis.longitudes(bias_maps[col])
        lats = ClimaAnalysis.latitudes(bias_maps[col])
        data = bias_maps[col].data

        GeoMakie.surface!(ax, lons, lats, data;
            colormap = cmap, colorrange = cmap_extrema,
            shading = GeoMakie.NoShading, rasterize = 2)
        GeoMakie.lines!(ax, GeoMakie.coastlines(); color = :black, linewidth = 0.5)

        # Add bias value and % change from 2010
        bias_val = biases[col]
        if col == 1
            label_text = @sprintf("Bias = %.3f", bias_val)
        else
            pct_change = 100 * (bias_val - biases[1]) / abs(biases[1])
            sign_str = pct_change >= 0 ? "+" : ""
            label_text = @sprintf("Bias = %.3f (%s%.1f%%)", bias_val, sign_str, pct_change)
        end
        GeoMakie.Label(fig[row, col, GeoMakie.Bottom()], label_text, fontsize = 20)
    end

    # Row label
    GeoMakie.Label(fig[row, 0],
        "$sn @ $(Int(PLOT_PRESSURE_LEVEL)) hPa",
        fontsize = 24,
        rotation = π/2)

    # Colorbar
    GeoMakie.Colorbar(fig[row, n_years + 1],
        colormap = cmap,
        colorrange = cmap_extrema,
        label = units_label,
        width = 15,
        ticklabelsize = 20,
        labelsize = 20,
        labelpadding = -5)

    @info "  $sn: Bias $(round(biases[1], sigdigits=3)) → $(round(biases[2], sigdigits=3)) → $(round(biases[3], sigdigits=3))"
end

# Set layout sizes
GeoMakie.colsize!(fig.layout, 0, GeoMakie.Fixed(70))   # Row labels
for col in 1:n_years
    GeoMakie.colsize!(fig.layout, col, GeoMakie.Fixed(370))
end
GeoMakie.colsize!(fig.layout, n_years + 1, GeoMakie.Fixed(65))  # Colorbar

GeoMakie.rowsize!(fig.layout, 0, GeoMakie.Fixed(55))   # Title row
for r in 1:n_vars
    GeoMakie.rowsize!(fig.layout, r, GeoMakie.Fixed(230))
end

GeoMakie.colgap!(fig.layout, -5)
GeoMakie.rowgap!(fig.layout, -5)

GeoMakie.resize_to_layout!(fig)

# Save
output_path_png = joinpath(plots_dir, "bias_comparison_clim_$(exp_suffix).png")
output_path_pdf = joinpath(plots_dir, "bias_comparison_clim_$(exp_suffix).pdf")
GeoMakie.save(output_path_png, fig)
GeoMakie.save(output_path_pdf, fig)
@info "Saved: $output_path_png"
@info "Saved: $output_path_pdf"

@info """
=====================================
Bias Comparison Complete (Climatology Reference)!
=====================================
Columns: $(join(years, ", "))
Reference: ERA5 February climatology (multi-year average)
Pressure level: $(Int(PLOT_PRESSURE_LEVEL)) hPa
Using: Ensemble means

Figures saved to:
  $output_path_png
  $output_path_pdf
"""
