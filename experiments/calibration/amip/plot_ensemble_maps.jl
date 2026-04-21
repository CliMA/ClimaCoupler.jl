#=
Plot tas and mslp maps from all ensemble members in iteration 0
to visualize the spatial patterns and ensemble spread.

Usage:
    julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/plot_ensemble_maps.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "ClimaEarth"))

using CairoMakie
using ClimaAnalysis
using GeoMakie
using Statistics
using JLD2
using Dates

# Load ClimaCoupler for CERESDataLoader
using ClimaCoupler
include(joinpath(pkgdir(ClimaCoupler), "src/CalibrationTools.jl"))

# Configuration
# const EXP_DIR = "/glade/derecho/scratch/cchristo/calibration/exp25"
const EXP_DIR = "/glade/derecho/scratch/zhaoyi/calibration/amip/exp1"
const EXP_NAME = basename(EXP_DIR)  # Extract exp name from path
const ITERATION = 2
const N_MEMBERS = 11
# Spinup and extension periods (must match run_calibration.jl)
const SPINUP = Dates.Day(7)
const EXTEND = Dates.Day(1)
# Date range for the calibration period (auto-inferred from model output if nothing)
const DATE_RANGE = (Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 31))  # e.g., (Dates.DateTime(2010, 1, 1), Dates.DateTime(2010, 1, 31))

# Fixed colorbar ranges for cross-iteration comparison
const TAS_MEAN_RANGE = (220.0, 310.0)    # K
const TAS_STD_RANGE = (0.0, 3.5)         # K
const MSLP_MEAN_RANGE = (960.0, 1050.0)  # hPa
const MSLP_STD_RANGE = (0.0, 3.5)        # hPa
const RSUT_MEAN_RANGE = (0.0, 350.0)     # W/m²
const RSUT_STD_RANGE = (0.0, 10.0)       # W/m²
const RLUT_MEAN_RANGE = (100.0, 320.0)   # W/m²
const RLUT_STD_RANGE = (0.0, 10.0)       # W/m²
const SWCRE_MEAN_RANGE = (-120.0, 0.0)   # W/m² (typically negative)
const SWCRE_STD_RANGE = (0.0, 10.0)      # W/m²
const LWCRE_MEAN_RANGE = (0.0, 80.0)     # W/m² (typically positive)
const LWCRE_STD_RANGE = (0.0, 10.0)      # W/m²

# Pressure level variable ranges
const TA_850_MEAN_RANGE = (240.0, 300.0)   # K at 850 hPa
const TA_500_MEAN_RANGE = (220.0, 270.0)   # K at 500 hPa
const TA_200_MEAN_RANGE = (200.0, 240.0)   # K at 200 hPa
const TA_STD_RANGE = (0.0, 5.0)            # K
const HUR_MEAN_RANGE = (0.0, 1.0)          # fractional (0-1)
const HUR_STD_RANGE = (0.0, 0.2)           # fractional

# Common pressure level variable names (for convenience)
const PRESSURE_LEVEL_VARS = [
    #"ta_850hPa", 
    "ta_500hPa", "ta_200hPa",
    #"hur_850hPa", "hur_500hPa", "hur_200hPa",
]

# Variables to plot (can include pressure level vars like "ta_850hPa")
PLOT_VARS = [
    #"ta_850hPa", 
    #"ta_500hPa", "ta_200hPa",
    #"hur_850hPa", "hur_500hPa", "hur_200hPa",
    "swcre", "lwcre",
]

# CERES variables (loaded directly from CERES, not from preprocessed_vars.jld2)
const CERES_VARS = Set(["rsut", "rlut", "rsutcs", "rlutcs", "rsds", "rsus", "rlds", "rlus"])

# Derived cloud radiative effect variables (computed from CERES components)
const CRE_VARS = Dict(
    "swcre" => ("rsutcs", "rsut"),   # SWCRE = rsutcs - rsut
    "lwcre" => ("rlutcs", "rlut"),   # LWCRE = rlutcs - rlut
)

"""
    is_cre_variable(short_name)

Check if a short_name is a derived cloud radiative effect variable.
"""
is_cre_variable(short_name::String) = haskey(CRE_VARS, short_name)

# Output directory for plots (created if it doesn't exist)
const OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "..", "plots", "ensemble_maps")

function ensure_output_dir()
    if !isdir(OUTPUT_DIR)
        mkpath(OUTPUT_DIR)
        println("Created output directory: $OUTPUT_DIR")
    end
    return OUTPUT_DIR
end

# ============================================================================
# Pressure Level Variable Support
# ============================================================================

"""
    is_pressure_level_variable(short_name)

Check if a short_name represents a pressure-level variable (e.g., "ta_850hPa").
"""
function is_pressure_level_variable(short_name::String)
    return occursin(r"_\d+hPa$", short_name)
end

"""
    parse_pressure_level_variable(short_name)

Parse a pressure-level variable name into (base_name, pressure_level).
e.g., "ta_850hPa" -> ("ta", 850.0)
"""
function parse_pressure_level_variable(short_name::String)
    m = match(r"^(.+)_(\d+)hPa$", short_name)
    if isnothing(m)
        error("Invalid pressure-level variable name: $short_name")
    end
    return (String(m.captures[1]), parse(Float64, m.captures[2]))
end

"""
    slice_pressure_level(var, pressure_hPa)

Slice an `OutputVar` at the given pressure level (in hPa), converting to Pa
if the data's pressure dimension is stored in Pa.

Uses ClimaAnalysis's built-in `pressure_name` to resolve the dimension
(handles both "pfull" and "pressure_level"), and `slice` already aliases
between these names, so we just need to pass the correct numeric value.
"""
function slice_pressure_level(var, pressure_hPa)
    pdim = ClimaAnalysis.pressure_name(var)
    dim_unit = get(get(var.dim_attributes, pdim, Dict()), "units", "")

    if dim_unit == "Pa"
        pressure_val = pressure_hPa * 100.0
    elseif dim_unit == "hPa"
        pressure_val = pressure_hPa
    else
        error("Unknown pressure units '$dim_unit' for dimension '$pdim'. " *
              "Expected 'Pa' or 'hPa'.")
    end

    return ClimaAnalysis.slice(var, pfull = pressure_val, by = ClimaAnalysis.MatchValue())
end

"""
Load raw ERA5 pressure level data for a variable.
Selects the month containing the calibration date range.
"""
function load_raw_era5_pressure_level(short_name; target_lons=nothing, target_lats=nothing, 
                                       date_range=nothing)
    base_name, pressure_hPa = parse_pressure_level_variable(short_name)
    
    loader = ClimaCoupler.CalibrationTools.ERA5PressureLevelDataLoader()
    var = Base.get(loader, base_name)
    
    # ERA5 artifact uses "pressure_level" dimension in Pa, convert hPa to Pa for slicing
    pressure_Pa = pressure_hPa * 100.0
    var = ClimaAnalysis.slice(var, pressure_level = pressure_Pa)
    
    # Select the month containing the calibration period (use end date to exclude spinup)
    if !isnothing(date_range)
        target_month = Dates.firstdayofmonth(date_range[2])
        println("  Selecting ERA5 $short_name for $(Dates.monthname(target_month)) $(Dates.year(target_month))")
        var = ClimaAnalysis.select(var, by = ClimaAnalysis.MatchValue(), time = target_month)
    end
    
    # Resample to target grid if provided
    if !isnothing(target_lons) && !isnothing(target_lats)
        var = ClimaAnalysis.resampled_as(var; lon=collect(target_lons), lat=collect(target_lats))
    end
    
    return var
end

# ============================================================================
# Date Range / Spinup Handling
# ============================================================================

"""
    infer_calibration_date_range(exp_dir, iteration, short_name)

Infer the calibration date range (excluding spinup and extension) from model output.

The simulation starts at `start_date` (which includes spinup before the actual
calibration period) and runs for `max(times)` seconds (which includes extension
after the calibration period). The calibration period is:

    [start_date + SPINUP,  start_date + max(times) - EXTEND]

This matches how model_interface.jl computes the simulation window:
    sim_start = calib_start - spinup
    sim_end   = calib_end   + extend
"""
function infer_calibration_date_range(exp_dir, iteration, short_name)
    # Load one member without windowing to read time metadata
    first_var = load_member_data(exp_dir, iteration, 1, short_name)
    times = ClimaAnalysis.times(first_var)

    if isempty(times)
        @warn "No time data found in model output"
        return nothing
    end

    start_date_attr = get(first_var.attributes, "start_date", nothing)
    if isnothing(start_date_attr)
        @warn "No start_date attribute found in model output"
        return nothing
    end

    sim_start = Dates.DateTime(start_date_attr)
    sim_end   = sim_start + Dates.Second(round(Int, maximum(times)))

    # Reverse the offsets applied in model_interface.jl
    calib_start = sim_start + SPINUP
    calib_end   = sim_end   - EXTEND

    println("  Simulation range:    $sim_start  to  $sim_end")
    println("  Calibration range:   $calib_start  to  $calib_end  (excluding $(SPINUP) spinup, $(EXTEND) extend)")

    return (calib_start, calib_end)
end

# ============================================================================
# Data Loading Functions
# ============================================================================

"""
Load raw (un-normalized) CERES data for a variable.
Selects the month containing the calibration date range.

CERES data is monthly with timestamps at the first of each month.
We use `select` with `MatchValue` to get exactly the month containing
the end of the date range (which excludes spinup).
"""
function load_raw_ceres(short_name; target_lons=nothing, target_lats=nothing, 
                        date_range=nothing)
    loader = ClimaCoupler.CalibrationTools.CERESDataLoader()
    var = Base.get(loader, short_name)
    
    # Select the month containing the calibration period (use end date to exclude spinup)
    if !isnothing(date_range)
        # CERES timestamps are at first of month, so select that month
        target_month = Dates.firstdayofmonth(date_range[2])
        println("  Selecting CERES $short_name for $(Dates.monthname(target_month)) $(Dates.year(target_month))")
        var = ClimaAnalysis.select(var, by = ClimaAnalysis.MatchValue(), time = target_month)
    end
    
    # Resample to target grid if provided
    if !isnothing(target_lons) && !isnothing(target_lats)
        var = ClimaAnalysis.resampled_as(var; lon=collect(target_lons), lat=collect(target_lats))
    end
    return var
end

"""
    load_member_data(exp_dir, iteration, member, short_name; date_range=nothing)

Load model output for a single ensemble member.
Handles pressure-level variables (e.g., "ta_850hPa") by loading from pressure coordinates.
If `date_range` is provided, window the output to exclude spinup time.
"""
function load_member_data(exp_dir, iteration, member, short_name; date_range=nothing)
    member_path = joinpath(exp_dir, "iteration_$(lpad(iteration, 3, '0'))", "member_$(lpad(member, 3, '0'))")
    
    # Find the job directory (e.g., wxquest_diagedmf)
    subdirs = filter(isdir, [joinpath(member_path, f) for f in readdir(member_path)])
    if isempty(subdirs)
        error("No output directory found in $member_path")
    end
    job_dir = first(subdirs)
    
    # SimDir expects the clima_atmos subdirectory
    atmos_path = joinpath(job_dir, "output_active", "clima_atmos")
    simdir = ClimaAnalysis.SimDir(atmos_path)
    
    # Handle pressure-level variables (e.g., "ta_850hPa")
    if is_pressure_level_variable(short_name)
        base_name, pressure_hPa = parse_pressure_level_variable(short_name)
        
        # Load from pressure coordinates and slice to level
        var = ClimaAnalysis.get(simdir; short_name=base_name, coord_type="pressure")
        var = slice_pressure_level(var, pressure_hPa)
        
        # Handle relative humidity units (convert % to fractional if needed)
        if base_name == "hur"
            units = ClimaAnalysis.units(var)
            if units == "%"
                var = ClimaAnalysis.convert_units(
                    var, "unitless"; 
                    conversion_function = x -> 0.01 * x
                )
            elseif units == "" || units == "1"
                var = ClimaAnalysis.set_units(var, "unitless")
            end
        end
    elseif is_cre_variable(short_name)
        # Cloud radiative effect: compute as clear-sky minus all-sky
        cs_name, as_name = CRE_VARS[short_name]
        var_cs = ClimaAnalysis.get(simdir; short_name=cs_name)
        var_as = ClimaAnalysis.get(simdir; short_name=as_name)
        var = var_cs - var_as
    else
        # Regular surface/2D variable
        var = ClimaAnalysis.get(simdir; short_name=short_name)
    end
    
    # Window to date_range to exclude spinup period if provided
    if !isnothing(date_range)
        var = ClimaAnalysis.window(var, "time"; left=date_range[1], right=date_range[2])
    end
    
    return var
end

"""
Get variable metadata for plotting (units, colormap, ranges, scale factor).
Handles both surface variables and pressure-level variables.
"""
function get_var_metadata(short_name)
    # Static metadata for known surface variables
    metadata = Dict(
        "tas" => (unit="K", cmap=:turbo, mean_range=TAS_MEAN_RANGE, std_range=TAS_STD_RANGE, scale=1.0, long_name="2m Temperature"),
        "mslp" => (unit="hPa", cmap=:turbo, mean_range=MSLP_MEAN_RANGE, std_range=MSLP_STD_RANGE, scale=0.01, long_name="Mean Sea Level Pressure"),
        "rsut" => (unit="W/m²", cmap=:YlOrRd, mean_range=RSUT_MEAN_RANGE, std_range=RSUT_STD_RANGE, scale=1.0, long_name="TOA Upward SW Radiation"),
        "rlut" => (unit="W/m²", cmap=:YlOrRd, mean_range=RLUT_MEAN_RANGE, std_range=RLUT_STD_RANGE, scale=1.0, long_name="TOA Upward LW Radiation"),
        # Pressure level variables
        "ta_850hPa" => (unit="K", cmap=:turbo, mean_range=TA_850_MEAN_RANGE, std_range=TA_STD_RANGE, scale=1.0, long_name="Temperature at 850 hPa"),
        "ta_500hPa" => (unit="K", cmap=:turbo, mean_range=TA_500_MEAN_RANGE, std_range=TA_STD_RANGE, scale=1.0, long_name="Temperature at 500 hPa"),
        "ta_200hPa" => (unit="K", cmap=:turbo, mean_range=TA_200_MEAN_RANGE, std_range=TA_STD_RANGE, scale=1.0, long_name="Temperature at 200 hPa"),
        "hur_850hPa" => (unit="", cmap=:Blues, mean_range=HUR_MEAN_RANGE, std_range=HUR_STD_RANGE, scale=1.0, long_name="Relative Humidity at 850 hPa"),
        "hur_500hPa" => (unit="", cmap=:Blues, mean_range=HUR_MEAN_RANGE, std_range=HUR_STD_RANGE, scale=1.0, long_name="Relative Humidity at 500 hPa"),
        "hur_200hPa" => (unit="", cmap=:Blues, mean_range=HUR_MEAN_RANGE, std_range=HUR_STD_RANGE, scale=1.0, long_name="Relative Humidity at 200 hPa"),
        # Derived CRE variables
        "swcre" => (unit="W/m²", cmap=:RdBu, mean_range=SWCRE_MEAN_RANGE, std_range=SWCRE_STD_RANGE, scale=1.0, long_name="SW Cloud Radiative Effect"),
        "lwcre" => (unit="W/m²", cmap=:RdBu, mean_range=LWCRE_MEAN_RANGE, std_range=LWCRE_STD_RANGE, scale=1.0, long_name="LW Cloud Radiative Effect"),
    )
    
    # Check for exact match first
    if haskey(metadata, short_name)
        return metadata[short_name]
    end
    
    # For unknown pressure-level variables, generate reasonable defaults
    if is_pressure_level_variable(short_name)
        base_name, pressure_hPa = parse_pressure_level_variable(short_name)
        if base_name == "ta"
            # Temperature defaults based on pressure level
            mean_range = pressure_hPa >= 800 ? TA_850_MEAN_RANGE : 
                         pressure_hPa >= 400 ? TA_500_MEAN_RANGE : TA_200_MEAN_RANGE
            return (unit="K", cmap=:turbo, mean_range=mean_range, std_range=TA_STD_RANGE, 
                    scale=1.0, long_name="Temperature at $(Int(pressure_hPa)) hPa")
        elseif base_name == "hur"
            return (unit="", cmap=:Blues, mean_range=HUR_MEAN_RANGE, std_range=HUR_STD_RANGE, 
                    scale=1.0, long_name="Relative Humidity at $(Int(pressure_hPa)) hPa")
        elseif base_name == "hus"
            return (unit="kg/kg", cmap=:YlGnBu, mean_range=(0.0, 0.02), std_range=(0.0, 0.005), 
                    scale=1.0, long_name="Specific Humidity at $(Int(pressure_hPa)) hPa")
        end
    end
    
    # Fallback for unknown variables
    return (unit="", cmap=:viridis, mean_range=(0.0, 1.0), std_range=(0.0, 1.0), scale=1.0, long_name=short_name)
end

"""
Time-average 3D data (assumes smallest dimension is time).
For model output with multiple daily snapshots, this computes the time mean.
For CERES monthly data, this just removes the singleton time dimension.
Replaces NaN values with 0 before averaging.
"""
function time_mean(data; verbose=false)
    # Replace NaN with 0 for averaging
    data_clean = replace(data, NaN => 0.0)
    
    if ndims(data_clean) == 3
        dims_sizes = size(data_clean)
        time_dim = argmin(dims_sizes)
        if verbose
            println("  Averaging over dim $time_dim (size $(dims_sizes[time_dim]))")
        end
        result = dropdims(mean(data_clean, dims=time_dim), dims=time_dim)
    else
        result = data_clean
    end
    
    # Replace any remaining NaN with 0
    return replace(result, NaN => 0.0)
end

"""
    plot_ensemble_maps(; exp_dir, iteration, n_members, short_names, date_range)

Plot ensemble mean and std maps for each variable.
If `date_range` is provided, window model output to exclude spinup period.
"""
function plot_ensemble_maps(; exp_dir=EXP_DIR, iteration=ITERATION, n_members=N_MEMBERS, 
                             short_names=PLOT_VARS, date_range=DATE_RANGE)
    println("Loading ensemble data from $exp_dir, iteration $iteration...")
    
    # Infer calibration date range (excluding spinup) if not provided
    if isnothing(date_range)
        try
            date_range = infer_calibration_date_range(exp_dir, iteration, short_names[1])
        catch e
            @warn "Could not infer calibration date range: $e"
        end
    end
    if !isnothing(date_range)
        println("  Windowing model output to: $(date_range[1]) - $(date_range[2])")
    end
    
    # Load all variables for all members
    all_vars = Dict{String, Vector{Any}}()
    for sn in short_names
        all_vars[sn] = []
    end
    
    n_loaded = 0
    for m in 1:n_members
        println("  Loading member $m...")
        try
            for sn in short_names
                var = load_member_data(exp_dir, iteration, m, sn; date_range=date_range)
                push!(all_vars[sn], var)
            end
            n_loaded += 1
        catch e
            println("  Warning: Failed to load member $m: $e")
        end
    end
    
    println("Loaded $n_loaded members successfully")
    
    if n_loaded == 0
        error("No members loaded!")
    end
    
    # Get lon/lat from first variable of first member
    first_var = all_vars[short_names[1]][1]
    lons = first_var.dims["lon"]
    lats = first_var.dims["lat"]
    
    # Debug info
    println("Data shape: $(size(first_var.data))")
    println("Dims available: $(keys(first_var.dims))")
    
    # Process each variable
    processed = Dict{String, NamedTuple}()
    for sn in short_names
        meta = get_var_metadata(sn)
        
        # Time-average all members (time_mean now handles NaN)
        all_data = [time_mean(v.data, verbose=(sn == short_names[1])) .* meta.scale for v in all_vars[sn]]
        
        # Filter out any members that are all zeros (failed runs)
        valid_data = filter(d -> any(d .!= 0.0), all_data)
        if isempty(valid_data)
            @warn "All members for $sn have zero/NaN data, using raw data anyway"
            valid_data = all_data
        elseif length(valid_data) < length(all_data)
            println("  $sn: Using $(length(valid_data))/$(length(all_data)) valid members")
        end
        
        # Compute statistics using valid members only
        data_min = minimum(minimum.(valid_data))
        data_max = maximum(maximum.(valid_data))
        data_mean = mean(valid_data)
        data_stack = cat(valid_data..., dims=3)
        data_std = dropdims(std(data_stack, dims=3), dims=3)
        
        # Replace any remaining NaN/Inf in statistics
        data_mean = replace(data_mean, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
        data_std = replace(data_std, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
        
        # Handle case where min/max are still NaN
        if isnan(data_min) || isnan(data_max)
            data_min = meta.mean_range[1]
            data_max = meta.mean_range[2]
            @warn "$sn: Using default range due to NaN values"
        end
        
        # Compute std range from data (use 2nd and 98th percentiles for robustness)
        std_flat = filter(isfinite, vec(data_std))
        if !isempty(std_flat)
            std_min = 0.0  # Always start at 0 for std
            std_max = quantile(std_flat, 0.98)  # 98th percentile to avoid outliers
            # Ensure reasonable minimum range
            if std_max < 0.1
                std_max = maximum(std_flat)
            end
        else
            std_min, std_max = meta.std_range
        end
        
        println("$sn range: $(round(data_min, digits=1)) to $(round(data_max, digits=1)) $(meta.unit), std max: $(round(std_max, digits=2))")
        
        processed[sn] = (all_data=valid_data, data_mean=data_mean, data_std=data_std, 
                         data_min=data_min, data_max=data_max, std_range=(std_min, std_max), meta=meta)
    end
    
    # Determine if we need to transpose for surface plot (expects lon x lat)
    needs_transpose = size(processed[short_names[1]].data_mean, 1) != length(lons)
    
    # Coastline data and colormaps
    coastline = GeoMakie.coastlines()
    std_cmap = cgrad([:white, :mediumpurple, :magenta, :orange, :yellow])
    
    # Create main figure with mean/std panels for each variable
    n_vars = length(short_names)
    fig = Figure(size=(1400, 500 * n_vars))
    
    for (row, sn) in enumerate(short_names)
        p = processed[sn]
        meta = p.meta
        
        plot_mean = needs_transpose ? permutedims(p.data_mean) : p.data_mean
        plot_std = needs_transpose ? permutedims(p.data_std) : p.data_std
        
        # Ensure colorranges are valid (not NaN)
        mean_range = meta.mean_range
        std_range = p.std_range  # Use data-derived std range
        if any(isnan.(mean_range))
            mean_range = (0.0, 1.0)
        end
        if any(isnan.(std_range)) || std_range[2] <= std_range[1]
            std_range = (0.0, 1.0)
        end
        
        # Mean panel
        ax1 = GeoAxis(fig[row, 1]; dest="+proj=robin",
            title="$(meta.long_name) Ensemble Mean ($(meta.unit)) - Iter $iteration")
        hm1 = surface!(ax1, lons, lats, plot_mean; colormap=meta.cmap, colorrange=mean_range, shading=NoShading)
        lines!(ax1, coastline, color=:black, linewidth=0.5)
        Colorbar(fig[row, 2], hm1)
        
        # Std panel
        ax2 = GeoAxis(fig[row, 3]; dest="+proj=robin",
            title="$(meta.long_name) Ensemble Std ($(meta.unit))")
        hm2 = surface!(ax2, lons, lats, plot_std; colormap=std_cmap, colorrange=std_range, shading=NoShading)
        lines!(ax2, coastline, color=:black, linewidth=0.5)
        Colorbar(fig[row, 4], hm2)
    end
    
    # Save main figure
    output_dir = ensure_output_dir()
    output_file = joinpath(output_dir, "ensemble_maps_$(EXP_NAME)_iter$(iteration).png")
    save(output_file, fig)
    println("Saved: $output_file")
    
    # Create individual member plots for each variable
    n_cols = 4
    n_rows = ceil(Int, n_loaded / n_cols)
    
    for sn in short_names
        p = processed[sn]
        meta = p.meta
        
        # Ensure valid colorrange for member plots
        member_min = isnan(p.data_min) ? meta.mean_range[1] : p.data_min
        member_max = isnan(p.data_max) ? meta.mean_range[2] : p.data_max
        if member_min >= member_max
            member_min, member_max = meta.mean_range
        end
        
        n_valid_members = length(p.all_data)
        n_rows_actual = ceil(Int, n_valid_members / n_cols)
        
        fig_members = Figure(size=(1600, 400 * max(n_rows_actual, 1)))
        Label(fig_members[0, 1:n_cols], "$(meta.long_name) for $n_valid_members valid ensemble members - Iteration $iteration", fontsize=18)
        
        for (i, data) in enumerate(p.all_data)
            row = div(i - 1, n_cols) + 1
            col = mod(i - 1, n_cols) + 1
            ax = GeoAxis(fig_members[row, col]; dest="+proj=robin",
                title="Member $i")
            plot_data = needs_transpose ? permutedims(data) : data
            # Replace any NaN in plot data
            plot_data = replace(plot_data, NaN => 0.0)
            surface!(ax, lons, lats, plot_data; colormap=meta.cmap, colorrange=(member_min, member_max), shading=NoShading)
            lines!(ax, coastline, color=:black, linewidth=0.3)
        end
        
        Colorbar(fig_members[1:max(n_rows_actual, 1), n_cols+1], limits=(member_min, member_max), colormap=meta.cmap, label="$(uppercase(sn)) ($(meta.unit))")
        
        output_file_members = joinpath(output_dir, "$(sn)_members_$(EXP_NAME)_iter$(iteration).png")
        save(output_file_members, fig_members)
        println("Saved: $output_file_members")
    end
    
    return fig
end

"""
    plot_bias_maps(; exp_dir, iteration, n_members, era5_vars_path, short_names, date_range)

Plot bias maps showing the difference between ensemble mean and ERA5 observations.
Creates a multi-panel figure with bias for each variable.

The `date_range` is used to select the appropriate observation period. If `nothing`,
it will be inferred from the model output timestamps.
"""
function plot_bias_maps(; 
    exp_dir=EXP_DIR, 
    iteration=ITERATION, 
    n_members=N_MEMBERS,
    era5_vars_path=joinpath(@__DIR__, "preprocessed_vars.jld2"),
    short_names=PLOT_VARS,
    date_range=DATE_RANGE
)
    println("\nLoading observation data...")
    
    # Classify variables by data source
    ceres_short_names = [sn for sn in short_names if sn in CERES_VARS]
    cre_short_names = [sn for sn in short_names if is_cre_variable(sn)]
    pressure_level_short_names = [sn for sn in short_names if is_pressure_level_variable(sn)]
    era5_surface_short_names = [sn for sn in short_names if !(sn in CERES_VARS) && !is_pressure_level_variable(sn) && !is_cre_variable(sn)]
    
    obs_data = Dict{String, Any}()
    
    # Infer calibration date range (excluding spinup) if not provided
    inferred_date_range = date_range
    if isnothing(inferred_date_range)
        try
            inferred_date_range = infer_calibration_date_range(exp_dir, iteration, short_names[1])
        catch e
            @warn "Could not infer calibration date range: $e"
        end
    end
    
    # Load ERA5 surface data from preprocessed file (if any non-CERES, non-pressure-level variables)
    if !isempty(era5_surface_short_names) && isfile(era5_vars_path)
        println("  Loading ERA5 surface variables from $era5_vars_path: $(join(era5_surface_short_names, ", "))")
        era5_vars = JLD2.load_object(era5_vars_path)
        for (key, var) in era5_vars
            sn = key[1]
            if sn in era5_surface_short_names && !haskey(obs_data, sn)
                # Check if data looks normalized (small range around 0)
                data_rng = maximum(var.data) - minimum(var.data)
                if data_rng < 20.0  # Likely normalized
                    @warn "  $sn data appears normalized (range=$data_rng), loading raw ERA5 instead"
                    # TODO: Load raw ERA5 for this variable
                else
                    obs_data[sn] = var
                end
            end
        end
    end
    
    # Load ERA5 pressure level data directly from artifact
    if !isempty(pressure_level_short_names)
        println("  Loading ERA5 pressure level variables: $(join(pressure_level_short_names, ", "))")
        for sn in pressure_level_short_names
            try
                var = load_raw_era5_pressure_level(sn; date_range=inferred_date_range)
                obs_data[sn] = var
                data_vals = filter(!isnan, vec(var.data))
                println("    Loaded $sn from ERA5 (range: $(round(minimum(data_vals), digits=2)) to $(round(maximum(data_vals), digits=2)))")
            catch e
                @warn "  Failed to load ERA5 pressure level $sn: $e"
            end
        end
    end
    
    # Load CERES data directly (always raw, never normalized)
    # NOTE: CERES data is monthly averages. If the model simulation is shorter than a month
    # (e.g., weekly runs), we're comparing sub-monthly model output to monthly CERES data.
    # This is an approximation - for rigorous comparisons, run the model for a full month.
    if !isempty(ceres_short_names)
        println("  Loading CERES variables directly: $(join(ceres_short_names, ", "))")
        if !isnothing(inferred_date_range)
            sim_days = Dates.value(inferred_date_range[2] - inferred_date_range[1]) / (1000 * 60 * 60 * 24)
            if sim_days < 28
                @warn "Model simulation is $(round(sim_days, digits=1)) days but CERES data is monthly. Bias comparison is approximate."
            end
        end
        for sn in ceres_short_names
            try
                var = load_raw_ceres(sn; date_range=inferred_date_range)
                obs_data[sn] = var
                println("    Loaded $sn from CERES (range: $(round(minimum(var.data), digits=1)) to $(round(maximum(var.data), digits=1)))")
            catch e
                @warn "  Failed to load CERES $sn: $e"
            end
        end
    end
    
    # Compute CRE variables from CERES components (clear-sky minus all-sky)
    if !isempty(cre_short_names)
        println("  Computing CRE variables from CERES: $(join(cre_short_names, ", "))")
        for sn in cre_short_names
            cs_name, as_name = CRE_VARS[sn]
            try
                var_cs = load_raw_ceres(cs_name; date_range=inferred_date_range)
                var_as = load_raw_ceres(as_name; date_range=inferred_date_range)
                # Resample to common grid if needed
                # if !isnothing(obs_lons) && !isnothing(obs_lats)
                #     var_cs = ClimaAnalysis.resampled_as(var_cs; lon=collect(obs_lons), lat=collect(obs_lats))
                #     var_as = ClimaAnalysis.resampled_as(var_as; lon=collect(obs_lons), lat=collect(obs_lats))
                # end
                obs_data[sn] = var_cs - var_as
                cre_vals = filter(!isnan, vec(obs_data[sn].data))
                println("    Computed $sn from CERES (range: $(round(minimum(cre_vals), digits=1)) to $(round(maximum(cre_vals), digits=1)))")
            catch e
                @warn "  Failed to compute $sn from CERES: $e"
            end
        end
    end
    
    # Check which variables we found
    found_vars = [sn for sn in short_names if haskey(obs_data, sn)]
    missing_vars = [sn for sn in short_names if !haskey(obs_data, sn)]
    
    if !isempty(missing_vars)
        println("  Warning: Missing observation data for: $(join(missing_vars, ", "))")
    end
    println("  Found observation data for: $(join(found_vars, ", "))")
    
    if isempty(found_vars)
        error("No observation variables found!")
    end
    
    # Get obs grid from first variable
    first_obs = obs_data[found_vars[1]]
    obs_lons = ClimaAnalysis.longitudes(first_obs)
    obs_lats = ClimaAnalysis.latitudes(first_obs)
    println("  Observation grid: $(length(obs_lons)) x $(length(obs_lats))")
    
    # Load ensemble data (windowed to exclude spinup if date_range provided)
    println("Loading ensemble data from $exp_dir, iteration $iteration...")
    if !isnothing(inferred_date_range)
        println("  Windowing model output to: $(inferred_date_range[1]) - $(inferred_date_range[2])")
    end
    model_vars = Dict{String, Vector{Any}}()
    for sn in found_vars
        model_vars[sn] = []
    end
    
    n_loaded = 0
    for m in 1:n_members
        try
            for sn in found_vars
                var = load_member_data(exp_dir, iteration, m, sn; date_range=inferred_date_range)
                push!(model_vars[sn], var)
            end
            n_loaded += 1
        catch e
            println("  Warning: Failed to load member $m: $e")
        end
    end
    
    println("  Loaded $n_loaded members")
    
    if n_loaded == 0
        error("No members loaded!")
    end
    
    # Get model lon/lat
    first_model = model_vars[found_vars[1]][1]
    model_lons = first_model.dims["lon"]
    model_lats = first_model.dims["lat"]
    
    # Resample observations to model grid if needed
    if length(model_lons) != length(obs_lons) || length(model_lats) != length(obs_lats)
        println("  Resampling observations to model grid ($(length(model_lons))x$(length(model_lats)))...")
        for sn in found_vars
            obs_data[sn] = ClimaAnalysis.resampled_as(obs_data[sn]; lon=collect(model_lons), lat=collect(model_lats))
        end
    end
    
    # Compute biases for each variable
    biases = Dict{String, NamedTuple}()
    needs_transpose = nothing
    
    for sn in found_vars
        meta = get_var_metadata(sn)
        
        # Get observation data with scaling (time-averaged to 2D)
        obs_raw = obs_data[sn].data
        obs_2d = time_mean(obs_raw)  # Remove time dim if present
        obs_scaled = obs_2d .* meta.scale
        
        # Compute ensemble mean (time-averaged) with scaling
        all_data = [time_mean(v.data) .* meta.scale for v in model_vars[sn]]
        model_mean = mean(all_data)
        
        # Check transpose on first variable
        if isnothing(needs_transpose)
            needs_transpose = size(model_mean, 1) != length(model_lons)
        end
        
        if needs_transpose
            model_mean = permutedims(model_mean)
        end
        
        # Ensure obs data has same orientation as model
        if needs_transpose && ndims(obs_scaled) == 2
            obs_scaled = permutedims(obs_scaled)
        end
        
        # Compute shared colorbar range from both obs and model data
        # Use 2nd and 98th percentiles for robustness to outliers
        all_values = vcat(vec(obs_scaled), vec(model_mean))
        valid_values = filter(isfinite, all_values)
        if !isempty(valid_values)
            shared_min = quantile(valid_values, 0.02)
            shared_max = quantile(valid_values, 0.98)
            # Ensure some minimum range
            if shared_max - shared_min < 1.0
                shared_min = minimum(valid_values)
                shared_max = maximum(valid_values)
            end
        else
            shared_min, shared_max = meta.mean_range
        end
        shared_range = (shared_min, shared_max)
        
        # Compute bias
        bias = model_mean .- obs_scaled
        bias_mean = mean(bias)
        bias_rmse = sqrt(mean(bias.^2))
        
        # Cap bias range for visibility
        bias_max = max(abs(minimum(bias)), abs(maximum(bias)))
        if sn in ["tas", "mslp"]
            bias_max = min(bias_max, 15.0)
        elseif sn in ["rsut", "rlut"]
            bias_max = min(bias_max, 50.0)  # W/m² cap for radiation
        elseif is_cre_variable(sn)
            bias_max = min(bias_max, 50.0)  # W/m² cap for CRE
        elseif is_pressure_level_variable(sn)
            base_name, _ = parse_pressure_level_variable(sn)
            if base_name == "ta"
                bias_max = min(bias_max, 8.0)  # K cap for pressure-level temperature
            elseif base_name == "hur"
                bias_max = min(bias_max, 0.2)   # fractional cap for relative humidity
            end
        end
        
        biases[sn] = (bias=bias, obs_data=obs_scaled, model_mean=model_mean, 
                      bias_mean=bias_mean, bias_rmse=bias_rmse, bias_max=bias_max, 
                      shared_range=shared_range, meta=meta)
        println("  $(uppercase(sn)) bias: mean=$(round(bias_mean, digits=2)) $(meta.unit), RMSE=$(round(bias_rmse, digits=2)) $(meta.unit)")
    end
    
    # Create figure with 2 columns, each variable gets 3 rows: Observations, Model, Bias
    n_vars = length(found_vars)
    n_cols = n_vars >= 2 ? 2 : 1  # Use 2 columns if we have 2+ variables
    n_var_rows = ceil(Int, n_vars / n_cols)  # Number of variable groups per column
    
    # Each variable needs 3 plot rows, figure width scales with columns
    # Column layout: [plot, colorbar, plot, colorbar] for 2 columns
    fig = Figure(size=(800 * n_cols, 350 * n_var_rows * 3))
    
    # Colormaps
    bias_cmap = Reverse(:RdBu)
    coastline = GeoMakie.coastlines()
    
    for (i, sn) in enumerate(found_vars)
        b = biases[sn]
        
        # Calculate position in 2-column grid
        var_idx = i - 1  # 0-indexed
        col_group = mod(var_idx, n_cols)  # Which column group (0 or 1)
        row_group = div(var_idx, n_cols)  # Which row group
        
        # Base column: each variable group takes 2 columns (plot + colorbar)
        base_col = col_group * 2 + 1
        # Base row: each variable takes 3 rows
        base_row = row_group * 3 + 1
        
        # Row 1: Observations (CERES/ERA5)
        ax_obs = GeoAxis(fig[base_row, base_col]; dest="+proj=robin",
            title="$(b.meta.long_name) - Observations (CERES/ERA5)")
        hm_obs = surface!(ax_obs, collect(model_lons), collect(model_lats), b.obs_data;
            colormap=b.meta.cmap, colorrange=b.shared_range, shading=NoShading)
        lines!(ax_obs, coastline, color=:black, linewidth=0.5)
        Colorbar(fig[base_row, base_col + 1], hm_obs, label=b.meta.unit)
        
        # Row 2: Model ensemble mean
        ax_model = GeoAxis(fig[base_row + 1, base_col]; dest="+proj=robin",
            title="$(b.meta.long_name) - Model Ensemble Mean")
        hm_model = surface!(ax_model, collect(model_lons), collect(model_lats), b.model_mean;
            colormap=b.meta.cmap, colorrange=b.shared_range, shading=NoShading)
        lines!(ax_model, coastline, color=:black, linewidth=0.5)
        Colorbar(fig[base_row + 1, base_col + 1], hm_model, label=b.meta.unit)
        
        # Row 3: Bias
        ax_bias = GeoAxis(fig[base_row + 2, base_col]; dest="+proj=robin",
            title="$(b.meta.long_name) Bias (Model - Obs) - Iter $iteration\nMean: $(round(b.bias_mean, digits=2)), RMSE: $(round(b.bias_rmse, digits=2)) $(b.meta.unit)")
        hm_bias = surface!(ax_bias, collect(model_lons), collect(model_lats), b.bias;
            colormap=bias_cmap, colorrange=(-b.bias_max, b.bias_max), shading=NoShading)
        lines!(ax_bias, coastline, color=:black, linewidth=0.5)
        Colorbar(fig[base_row + 2, base_col + 1], hm_bias, label=b.meta.unit)
    end
    
    # Save to dedicated plots directory
    output_dir = ensure_output_dir()
    output_file = joinpath(output_dir, "bias_$(EXP_NAME)_iter$(iteration).png")
    save(output_file, fig)
    println("Saved: $output_file")
    
    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    plot_ensemble_maps()
    plot_bias_maps()
end
