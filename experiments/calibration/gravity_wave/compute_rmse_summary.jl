#=
compute_rmse_summary.jl

Compute RMSE for each iteration against year-specific ERA5 data.

- Iteration 0 ensemble vs ERA5 2010 (initial params)
- Iteration 1 ensemble vs ERA5 2011 (after 1 calibration step)
- Iteration 2 ensemble vs ERA5 2012 (after 2 calibration steps)

Usage:
    julia --project=experiments/ClimaEarth compute_rmse_summary.jl [output_dir]
=#

using ClimaAnalysis
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

@info "RMSE Summary (year-specific ERA5 comparison)"
@info "Output directory: $output_dir"

short_names = ["ta", "ua", "va"]
iterations = [0, 1, 2]
years = [2010, 2011, 2012]

# ERA5 pressure levels and target level for RMSE
const ERA5_PRESSURE_LEVELS = [850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 250.0, 200.0]
const TARGET_PRESSURE = 500.0  # hPa

# Units for each variable
var_units = Dict(
    "ta" => "K",
    "ua" => "m s^-1",
    "va" => "m s^-1",
)

# ==========================================================================
# DATA LOADING FUNCTIONS
# ==========================================================================

"""
Load ERA5 observation variable and window to a specific year/month.
"""
function load_era5_var_for_date(short_name, year, month=2)
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

    # Window to specific date
    # ERA5 uses 1979-01-01 as epoch, time is in seconds
    target_date = Dates.DateTime(year, month, 1)
    era5_epoch = Dates.DateTime(1979, 1, 1)
    target_time = Float64(Dates.value(target_date - era5_epoch) / 1000)  # ms to seconds

    if ClimaAnalysis.has_time(era5_var)
        time_vals = ClimaAnalysis.times(era5_var)
        # Find closest time
        idx = argmin(abs.(time_vals .- target_time))
        era5_var = ClimaAnalysis.slice(era5_var; time = time_vals[idx])
    end

    return era5_var
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

    return ensemble_mean, length(member_vars)
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
# COMPUTE RMSE FOR EACH ITERATION vs YEAR-SPECIFIC ERA5
# ==========================================================================

println("\n" * "="^70)
println("RMSE Summary: Each iteration vs year-specific ERA5 @ $(Int(TARGET_PRESSURE)) hPa")
println("="^70)

# Store results
results = Dict{String, Vector{Float64}}()

for short_name in short_names
    results[short_name] = Float64[]

    for (i, iteration) in enumerate(iterations)
        year = years[i]

        # Load year-specific ERA5
        era5_var = load_era5_var_for_date(short_name, year, 2)  # February

        # Load simulation ensemble mean
        sim_var, n_members = load_ensemble_mean(output_dir, iteration, short_name)

        # Compute RMSE
        rmse = compute_rmse_at_level(sim_var, era5_var, TARGET_PRESSURE)
        push!(results[short_name], rmse)

        if short_name == short_names[1]
            @info "Iteration $iteration: simulated $(year), compared to ERA5 Feb $(year), n_members=$n_members"
        end
    end
end

# Print results table
println("\n" * "-"^70)
@printf("%-8s  %-20s  %-20s  %-20s\n",
    "Var",
    "Iter 0 (2010→ERA5'10)",
    "Iter 1 (2011→ERA5'11)",
    "Iter 2 (2012→ERA5'12)")
println("-"^70)

for short_name in short_names
    r = results[short_name]
    @printf("%-8s  %18.4f  %18.4f  %18.4f\n", short_name, r[1], r[2], r[3])
end

println("-"^70)

# Compute and print changes relative to iteration 0
println("\nChange relative to Iteration 0 (initial params):")
println("-"^70)
@printf("%-8s  %-20s  %-20s  %-20s\n", "Var", "Iter 0 (baseline)", "Iter 1 change", "Iter 2 change")
println("-"^70)

for short_name in short_names
    r = results[short_name]
    baseline = r[1]
    change1 = 100 * (r[2] - baseline) / baseline
    change2 = 100 * (r[3] - baseline) / baseline
    sign1 = change1 >= 0 ? "+" : ""
    sign2 = change2 >= 0 ? "+" : ""
    @printf("%-8s  %18.4f  %17s%.1f%%  %17s%.1f%%\n",
        short_name, baseline, sign1, change1, sign2, change2)
end

println("="^70)

println("\nNOTE: Each iteration simulates a DIFFERENT year:")
println("  - Iter 0: Initial params, simulates Feb 2010, compared to ERA5 Feb 2010")
println("  - Iter 1: Calibrated (1 step), simulates Feb 2011, compared to ERA5 Feb 2011")
println("  - Iter 2: Calibrated (2 steps), simulates Feb 2012, compared to ERA5 Feb 2012")
println("\nChanges reflect BOTH calibration effects AND interannual variability.")
println("To isolate calibration effect, would need to run calibrated params on 2010.")
