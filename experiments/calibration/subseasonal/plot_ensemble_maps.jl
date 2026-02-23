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

# Configuration
# const EXP_DIR = "/glade/derecho/scratch/cchristo/calibration/exp25"
const EXP_DIR = "/glade/derecho/scratch/zhaoyi/calibration/exp7"
const EXP_NAME = basename(EXP_DIR)  # Extract exp name from path
const ITERATION = 3
const N_MEMBERS = 1

# Fixed colorbar ranges for cross-iteration comparison
const TAS_MEAN_RANGE = (220.0, 310.0)    # K
const TAS_STD_RANGE = (0.0, 3.5)         # K
const MSLP_MEAN_RANGE = (960.0, 1050.0)  # hPa
const MSLP_STD_RANGE = (0.0, 3.5)        # hPa
const RSUT_MEAN_RANGE = (0.0, 350.0)     # W/m²
const RSUT_STD_RANGE = (0.0, 10.0)       # W/m²
const RLUT_MEAN_RANGE = (100.0, 320.0)   # W/m²
const RLUT_STD_RANGE = (0.0, 10.0)       # W/m²

# Variables to plot
const PLOT_VARS = ["rsut", "rlut"]

# Output directory for plots (created if it doesn't exist)
const OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "..", "plots", "ensemble_maps")

function ensure_output_dir()
    if !isdir(OUTPUT_DIR)
        mkpath(OUTPUT_DIR)
        println("Created output directory: $OUTPUT_DIR")
    end
    return OUTPUT_DIR
end

function load_member_data(exp_dir, iteration, member, short_name)
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
    
    # Get the variable
    var = ClimaAnalysis.get(simdir; short_name=short_name)
    return var
end

"""
Get variable metadata for plotting (units, colormap, ranges, scale factor).
"""
function get_var_metadata(short_name)
    metadata = Dict(
        "tas" => (unit="K", cmap=:turbo, mean_range=TAS_MEAN_RANGE, std_range=TAS_STD_RANGE, scale=1.0, long_name="2m Temperature"),
        "mslp" => (unit="hPa", cmap=:turbo, mean_range=MSLP_MEAN_RANGE, std_range=MSLP_STD_RANGE, scale=0.01, long_name="Mean Sea Level Pressure"),
        "rsut" => (unit="W/m²", cmap=:YlOrRd, mean_range=RSUT_MEAN_RANGE, std_range=RSUT_STD_RANGE, scale=1.0, long_name="TOA Upward SW Radiation"),
        "rlut" => (unit="W/m²", cmap=:YlOrRd, mean_range=RLUT_MEAN_RANGE, std_range=RLUT_STD_RANGE, scale=1.0, long_name="TOA Upward LW Radiation"),
    )
    return get(metadata, short_name, (unit="", cmap=:viridis, mean_range=(0.0, 1.0), std_range=(0.0, 1.0), scale=1.0, long_name=short_name))
end

"""
Time-average 3D data (assumes smallest dimension is time).
"""
function time_mean(data; verbose=false)
    if ndims(data) == 3
        dims_sizes = size(data)
        time_dim = argmin(dims_sizes)
        if verbose
            println("  Averaging over dim $time_dim (size $(dims_sizes[time_dim]))")
        end
        return dropdims(mean(data, dims=time_dim), dims=time_dim)
    else
        return data
    end
end

function plot_ensemble_maps(; exp_dir=EXP_DIR, iteration=ITERATION, n_members=N_MEMBERS, short_names=PLOT_VARS)
    println("Loading ensemble data from $exp_dir, iteration $iteration...")
    
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
                var = load_member_data(exp_dir, iteration, m, sn)
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
        
        # Time-average all members
        all_data = [time_mean(v.data, verbose=(sn == short_names[1])) .* meta.scale for v in all_vars[sn]]
        
        # Compute statistics
        data_min = minimum(minimum.(all_data))
        data_max = maximum(maximum.(all_data))
        data_mean = mean(all_data)
        data_stack = cat(all_data..., dims=3)
        data_std = dropdims(std(data_stack, dims=3), dims=3)
        
        println("$sn range: $(round(data_min, digits=1)) to $(round(data_max, digits=1)) $(meta.unit)")
        
        processed[sn] = (all_data=all_data, data_mean=data_mean, data_std=data_std, 
                         data_min=data_min, data_max=data_max, meta=meta)
    end
    
    # Determine if we need to transpose for heatmap (expects lon x lat)
    needs_transpose = size(processed[short_names[1]].data_mean, 1) != length(lons)
    
    # Get coastline data
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
        
        # Mean panel
        ax1 = Axis(fig[row, 1], 
            title="$(meta.long_name) Ensemble Mean ($(meta.unit)) - Iter $iteration",
            xlabel="Longitude", ylabel="Latitude")
        hm1 = heatmap!(ax1, lons, lats, plot_mean, colormap=meta.cmap, colorrange=meta.mean_range)
        lines!(ax1, coastline, color=:black, linewidth=0.5)
        Colorbar(fig[row, 2], hm1)
        
        # Std panel
        ax2 = Axis(fig[row, 3], 
            title="$(meta.long_name) Ensemble Std ($(meta.unit))",
            xlabel="Longitude", ylabel="Latitude")
        hm2 = heatmap!(ax2, lons, lats, plot_std, colormap=std_cmap, colorrange=meta.std_range)
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
        
        fig_members = Figure(size=(1600, 400 * n_rows))
        Label(fig_members[0, 1:n_cols], "$(meta.long_name) for all $n_loaded ensemble members - Iteration $iteration", fontsize=18)
        
        for (i, data) in enumerate(p.all_data)
            row = div(i - 1, n_cols) + 1
            col = mod(i - 1, n_cols) + 1
            ax = Axis(fig_members[row, col], title="Member $i", aspect=DataAspect())
            plot_data = needs_transpose ? permutedims(data) : data
            heatmap!(ax, lons, lats, plot_data, colormap=meta.cmap, colorrange=(p.data_min, p.data_max))
            lines!(ax, coastline, color=:black, linewidth=0.3)
            hidedecorations!(ax)
        end
        
        Colorbar(fig_members[1:n_rows, n_cols+1], limits=(p.data_min, p.data_max), colormap=meta.cmap, label="$(uppercase(sn)) ($(meta.unit))")
        
        output_file_members = joinpath(output_dir, "$(sn)_members_$(EXP_NAME)_iter$(iteration).png")
        save(output_file_members, fig_members)
        println("Saved: $output_file_members")
    end
    
    return fig
end

"""
    plot_bias_maps(; exp_dir, iteration, n_members, era5_vars_path, short_names)

Plot bias maps showing the difference between ensemble mean and ERA5 observations.
Creates a multi-panel figure with bias for each variable.
"""
function plot_bias_maps(; 
    exp_dir=EXP_DIR, 
    iteration=ITERATION, 
    n_members=N_MEMBERS,
    era5_vars_path=joinpath(@__DIR__, "preprocessed_vars.jld2"),
    short_names=PLOT_VARS
)
    println("\nLoading ERA5 observations from $era5_vars_path...")
    
    # Load preprocessed ERA5 variables
    era5_vars = JLD2.load_object(era5_vars_path)
    
    # Find requested variables in ERA5 data
    # Keys are (short_name, (start_date, end_date)) tuples
    era5_data = Dict{String, Any}()
    for (key, var) in era5_vars
        sn = key[1]
        if sn in short_names && !haskey(era5_data, sn)
            era5_data[sn] = var
        end
    end
    
    # Check which variables we found
    found_vars = [sn for sn in short_names if haskey(era5_data, sn)]
    missing_vars = [sn for sn in short_names if !haskey(era5_data, sn)]
    
    if !isempty(missing_vars)
        println("  Warning: Missing ERA5 data for: $(join(missing_vars, ", "))")
    end
    println("  Found ERA5 data for: $(join(found_vars, ", "))")
    
    if isempty(found_vars)
        error("No ERA5 variables found!")
    end
    
    # Get ERA5 grid from first variable
    first_era5 = era5_data[found_vars[1]]
    era5_lons = ClimaAnalysis.longitudes(first_era5)
    era5_lats = ClimaAnalysis.latitudes(first_era5)
    println("  ERA5 grid: $(length(era5_lons)) x $(length(era5_lats))")
    
    # Load ensemble data
    println("Loading ensemble data from $exp_dir, iteration $iteration...")
    model_vars = Dict{String, Vector{Any}}()
    for sn in found_vars
        model_vars[sn] = []
    end
    
    n_loaded = 0
    for m in 1:n_members
        try
            for sn in found_vars
                var = load_member_data(exp_dir, iteration, m, sn)
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
    
    # Resample ERA5 to model grid if needed
    if length(model_lons) != length(era5_lons) || length(model_lats) != length(era5_lats)
        println("  Resampling ERA5 to model grid ($(length(model_lons))x$(length(model_lats)))...")
        for sn in found_vars
            era5_data[sn] = ClimaAnalysis.resampled_as(era5_data[sn]; lon=collect(model_lons), lat=collect(model_lats))
        end
    end
    
    # Compute biases for each variable
    biases = Dict{String, NamedTuple}()
    needs_transpose = nothing
    
    for sn in found_vars
        meta = get_var_metadata(sn)
        
        # Get ERA5 data with scaling
        era5_scaled = era5_data[sn].data .* meta.scale
        
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
        
        # Compute bias
        bias = model_mean .- era5_scaled
        bias_mean = mean(bias)
        bias_rmse = sqrt(mean(bias.^2))
        
        # Cap bias range for visibility
        bias_max = max(abs(minimum(bias)), abs(maximum(bias)))
        if sn in ["tas", "mslp"]
            bias_max = min(bias_max, 15.0)
        elseif sn in ["rsut", "rlut"]
            bias_max = min(bias_max, 50.0)  # W/m² cap for radiation
        end
        
        biases[sn] = (bias=bias, bias_mean=bias_mean, bias_rmse=bias_rmse, bias_max=bias_max, meta=meta)
        println("  $(uppercase(sn)) bias: mean=$(round(bias_mean, digits=2)) $(meta.unit), RMSE=$(round(bias_rmse, digits=2)) $(meta.unit)")
    end
    
    # Create bias figure - 2 columns layout
    n_vars = length(found_vars)
    n_rows = ceil(Int, n_vars / 2)
    fig = Figure(size=(1400, 500 * n_rows))
    
    # Diverging colormap for bias
    bias_cmap = Reverse(:RdBu)
    coastline = GeoMakie.coastlines()
    
    for (i, sn) in enumerate(found_vars)
        b = biases[sn]
        row = div(i - 1, 2) + 1
        col_base = mod(i - 1, 2) * 2 + 1  # 1 or 3
        
        ax = Axis(fig[row, col_base], 
            title="$(b.meta.long_name) Bias (Model - ERA5) - Iter $iteration\nMean: $(round(b.bias_mean, digits=2)), RMSE: $(round(b.bias_rmse, digits=2)) $(b.meta.unit)",
            xlabel="Longitude", ylabel="Latitude")
        hm = heatmap!(ax, collect(model_lons), collect(model_lats), b.bias, 
            colormap=bias_cmap, colorrange=(-b.bias_max, b.bias_max))
        lines!(ax, coastline, color=:black, linewidth=0.5)
        Colorbar(fig[row, col_base + 1], hm, label=b.meta.unit)
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
