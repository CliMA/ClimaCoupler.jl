#=
plot_rmse_barchart.jl

Create a per-year RMSE comparison bar chart for gravity wave calibration.
Shows RMSE for each year (2010, 2011, 2012) grouped by variable (ta, ua, va).

This frames the calibration as sequential validation:
- Year 2010: Initial parameters
- Year 2011: After 1 calibration iteration
- Year 2012: After 2 calibration iterations

Usage:
    julia --project=experiments/ClimaEarth plot_rmse_barchart.jl [output_dir]
=#

using ClimaAnalysis
using CairoMakie
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

plots_dir = joinpath(output_dir, "plots")
mkpath(plots_dir)

@info "Per-Year RMSE Bar Chart (Gravity Wave Calibration)"
@info "Each iteration compared to ERA5 February climatology (multi-year average)"
@info "Output directory: $output_dir"

# Variables and years
short_names = ["ta", "ua", "va"]
years = [2010, 2011, 2012]
iterations = [0, 1, 2]  # iteration i corresponds to year[i+1]

# ERA5 pressure levels
const ERA5_PRESSURE_LEVELS = [850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 250.0, 200.0]

# Units for each variable
var_units = Dict(
    "ta" => "K",
    "ua" => "m s^-1",
    "va" => "m s^-1",
)

# Variable display names
var_names = Dict(
    "ta" => "Temperature",
    "ua" => "Zonal Wind",
    "va" => "Meridional Wind",
)

# ==========================================================================
# DATA LOADING FUNCTIONS (adapted from plot_rmse_comparison_gw.jl)
# ==========================================================================

# ERA5 epoch for time conversion
const ERA5_EPOCH = Dates.DateTime(1979, 1, 1)

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

    # Find all time points that correspond to the target month (e.g., February)
    time_vals = ClimaAnalysis.times(era5_var)
    feb_indices = Int[]

    for (i, t) in enumerate(time_vals)
        # Convert time (seconds since ERA5 epoch) to DateTime
        dt = ERA5_EPOCH + Dates.Second(round(Int, t))
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

    return climatology
end

"""
Load simulation variable and convert z-levels to pressure coordinates.
"""
function get_var(short_name, simdir)
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
function load_ensemble_mean(output_dir, iteration, short_name)
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

    return ensemble_mean
end

"""
Compute global RMSE (latitude-weighted) averaged over all pressure levels.
"""
function compute_global_rmse(sim_var, obs_var)
    diff = sim_var - obs_var
    sq_diff = ClimaAnalysis.replace(x -> x^2, diff)
    mean_sq = ClimaAnalysis.average_lonlat(sq_diff; weighted=true)

    # Average over pressure levels if present
    if ClimaAnalysis.has_pressure(mean_sq)
        mean_sq = ClimaAnalysis.average_lat(mean_sq)  # This averages over remaining dims
        return sqrt(mean(mean_sq.data))
    else
        return sqrt(mean_sq.data[1])
    end
end

"""
Compute global RMSE at a specific pressure level.
"""
function compute_rmse_at_level(sim_var, obs_var, pressure_level)
    sim_slice = ClimaAnalysis.slice(sim_var; pressure_level)
    obs_slice = ClimaAnalysis.slice(obs_var; pressure_level)

    # Resample to same grid
    obs_lons = ClimaAnalysis.longitudes(obs_slice)
    obs_lats = ClimaAnalysis.latitudes(obs_slice)
    sim_resampled = ClimaAnalysis.resampled_as(sim_slice; longitude=obs_lons, latitude=obs_lats)

    # Match units
    obs_units = ClimaAnalysis.units(obs_slice)
    sim_resampled = ClimaAnalysis.set_units(sim_resampled, obs_units)

    diff = sim_resampled - obs_slice
    sq_diff = ClimaAnalysis.replace(x -> x^2, diff)
    mean_sq = ClimaAnalysis.average_lonlat(sq_diff; weighted=true)
    return sqrt(mean_sq.data[1])
end

# ==========================================================================
# COMPUTE RMSE FOR ALL YEARS AND VARIABLES
# ==========================================================================

@info "Computing RMSE for each iteration vs ERA5 February climatology..."

# Store results: rmse_data[var][iteration] = rmse_value
rmse_data = Dict{String, Vector{Float64}}()

for short_name in short_names
    rmse_data[short_name] = Float64[]

    # Load ERA5 February climatology (same reference for all iterations)
    @info "Loading ERA5 climatology for $short_name..."
    era5_clim = load_era5_climatology(short_name, 2)  # February climatology

    for (i, iteration) in enumerate(iterations)
        year = years[i]
        @info "  Processing $short_name, iteration $iteration (simulated $year)..."

        # Load simulation ensemble mean
        sim_var = load_ensemble_mean(output_dir, iteration, short_name)

        # Compute RMSE at 500 hPa against climatology
        rmse = compute_rmse_at_level(sim_var, era5_clim, 500.0)
        push!(rmse_data[short_name], rmse)

        @info "    RMSE @ 500 hPa: $(round(rmse, sigdigits=3)) [vs Feb climatology]"
    end
end

# ==========================================================================
# CREATE BAR CHART WITH BASELINE REFERENCE
# ==========================================================================

@info "Creating bar chart with baseline reference..."

fig = Figure(size = (600, 210), fontsize = 12)
ax = Axis(fig[1, 1],
    ylabel = "RMSE",
    title = "RMSE vs ERA5 46-year February Climatology (500 hPa)",
    xticks = (1:3, ["2010", "2011", "2012"]),
    ylabelsize = 12,
    titlesize = 12,
)

# Bar positions and width
n_vars = length(short_names)
bar_width = 0.25
offsets = [-bar_width, 0, bar_width]

# Colors for each variable
colors = [:steelblue, :coral, :seagreen]

# Find max RMSE to set ylims with room for labels
max_rmse = maximum(maximum(rmse_data[sn]) for sn in short_names)
ylims!(ax, 0, max_rmse * 1.30)  # 30% headroom for top labels

# Plot bars for each variable
for (j, short_name) in enumerate(short_names)
    xs = (1:3) .+ offsets[j]
    ys = rmse_data[short_name]
    barplot!(ax, xs, ys,
        width = bar_width,
        color = colors[j],
        label = "$(var_names[short_name]) ($(short_name))"
    )

    # Add horizontal dashed line at 2010 baseline level
    baseline = rmse_data[short_name][1]
    hlines!(ax, [baseline], color = colors[j], linestyle = :dash, linewidth = 1.5, alpha = 0.6)

    # Add value labels on top of bars
    for (i, (x, y)) in enumerate(zip(xs, ys))
        pct_change = 100 * (y - baseline) / baseline
        if i == 1
            # Baseline year - just show value
            text!(ax, x, y + 0.15, text = @sprintf("%.2f", y),
                align = (:center, :bottom), fontsize = 10, font = :bold)
        else
            # Validation years - show value and % change
            sign_str = pct_change >= 0 ? "+" : ""
            text!(ax, x, y + 0.15, text = @sprintf("%.2f\n(%s%.1f%%)", y, sign_str, pct_change),
                align = (:center, :bottom), fontsize = 10)
        end
    end
end

# Add legend with annotation text
Legend(fig[1, 2], ax, "Variables\n\nDashed lines =\n2010 baseline", framevisible = true, labelsize = 10, titlesize = 10)

# Save PNG and PDF
output_path_png = joinpath(plots_dir, "rmse_per_year.png")
output_path_pdf = joinpath(plots_dir, "rmse_per_year.pdf")
save(output_path_png, fig)
save(output_path_pdf, fig)
@info "Saved: $output_path_png"
@info "Saved: $output_path_pdf"

# Print summary table
println("\n" * "="^70)
println("RMSE Summary (500 hPa) - vs ERA5 February Climatology")
println("="^70)
@printf("%-12s %18s %18s %18s\n", "Variable", "Iter 0 (Initial)", "Iter 1 (Calib)", "Iter 2 (Calib)")
println("-"^70)
for short_name in short_names
    r = rmse_data[short_name]
    baseline = r[1]
    change1 = 100 * (r[2] - baseline) / baseline
    change2 = 100 * (r[3] - baseline) / baseline
    sign1 = change1 >= 0 ? "+" : ""
    sign2 = change2 >= 0 ? "+" : ""
    @printf("%-12s %18.3f %12.3f (%s%.1f%%) %12.3f (%s%.1f%%)\n",
        short_name, r[1], r[2], sign1, change1, r[3], sign2, change2)
end
println("="^70)
println("\nNote: All iterations compared to the SAME ERA5 February climatology.")
println("  % change relative to Iter 0 (initial params) is meaningful.")
println("  - Iter 0: Initial params, simulated Feb 2010")
println("  - Iter 1: Calibrated params, simulated Feb 2011")
println("  - Iter 2: Further calibrated, simulated Feb 2012")
