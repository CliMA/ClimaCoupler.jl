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
using Statistics

# Configuration
const EXP_DIR = "/glade/derecho/scratch/cchristo/calibration/exp24"
const ITERATION = 0
const N_MEMBERS = 11

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

function plot_ensemble_maps(; exp_dir=EXP_DIR, iteration=ITERATION, n_members=N_MEMBERS)
    println("Loading ensemble data from $exp_dir, iteration $iteration...")
    
    # Load tas and mslp for all members
    tas_vars = []
    mslp_vars = []
    
    for m in 1:n_members
        println("  Loading member $m...")
        try
            tas = load_member_data(exp_dir, iteration, m, "tas")
            mslp = load_member_data(exp_dir, iteration, m, "mslp")
            push!(tas_vars, tas)
            push!(mslp_vars, mslp)
        catch e
            println("  Warning: Failed to load member $m: $e")
        end
    end
    
    n_loaded = length(tas_vars)
    println("Loaded $n_loaded members successfully")
    
    if n_loaded == 0
        error("No members loaded!")
    end
    
    # Get lon/lat from first member
    lons = tas_vars[1].dims["lon"]
    lats = tas_vars[1].dims["lat"]
    
    # Debug: print data shape
    println("Data shape: $(size(tas_vars[1].data))")
    println("Dims available: $(keys(tas_vars[1].dims))")
    
    # Average over time dimension if present
    # ClimaAnalysis data is typically (lon, lat, time) or (time, lon, lat)
    function time_mean(data)
        if ndims(data) == 3
            # Check which dimension is time (smallest one, likely 7 for weekly)
            dims_sizes = size(data)
            time_dim = argmin(dims_sizes)
            println("  Averaging over dim $time_dim (size $(dims_sizes[time_dim]))")
            return dropdims(mean(data, dims=time_dim), dims=time_dim)
        else
            return data
        end
    end
    
    # Compute global min/max for consistent colorbars
    all_tas = [time_mean(v.data) for v in tas_vars]
    all_mslp = [time_mean(v.data) for v in mslp_vars]
    
    tas_min = minimum(minimum.(all_tas))
    tas_max = maximum(maximum.(all_tas))
    mslp_min = minimum(minimum.(all_mslp)) / 100  # Convert to hPa
    mslp_max = maximum(maximum.(all_mslp)) / 100
    
    println("tas range: $(round(tas_min, digits=1)) to $(round(tas_max, digits=1)) K")
    println("mslp range: $(round(mslp_min, digits=1)) to $(round(mslp_max, digits=1)) hPa")
    
    # Compute ensemble mean (using time-averaged data)
    tas_mean = mean(all_tas)
    mslp_mean = mean(all_mslp) ./ 100  # hPa
    
    # Compute ensemble std (spread)
    tas_stack = cat(all_tas..., dims=3)
    mslp_stack = cat(all_mslp..., dims=3)
    tas_std = dropdims(std(tas_stack, dims=3), dims=3)
    mslp_std = dropdims(std(mslp_stack, dims=3), dims=3) ./ 100  # hPa
    
    # Create figure with 4 panels: tas mean, tas std, mslp mean, mslp std
    fig = Figure(size=(1400, 1000))
    
    # Determine if we need to transpose for heatmap (expects lon x lat)
    needs_transpose = size(tas_mean, 1) != length(lons)
    plot_tas_mean = needs_transpose ? permutedims(tas_mean) : tas_mean
    plot_tas_std = needs_transpose ? permutedims(tas_std) : tas_std
    plot_mslp_mean = needs_transpose ? permutedims(mslp_mean) : mslp_mean
    plot_mslp_std = needs_transpose ? permutedims(mslp_std) : mslp_std
    
    # Use turbo colormap (similar to gist_ncar from matplotlib)
    tas_cmap = :turbo
    
    # TAS mean
    ax1 = Axis(fig[1, 1], 
        title="TAS Ensemble Mean (K) - Iteration $iteration",
        xlabel="Longitude", ylabel="Latitude")
    hm1 = heatmap!(ax1, lons, lats, plot_tas_mean, colormap=tas_cmap)
    Colorbar(fig[1, 2], hm1)
    
    # TAS std
    ax2 = Axis(fig[1, 3], 
        title="TAS Ensemble Std (K)",
        xlabel="Longitude", ylabel="Latitude")
    hm2 = heatmap!(ax2, lons, lats, plot_tas_std, colormap=:viridis)
    Colorbar(fig[1, 4], hm2)
    
    # MSLP mean
    ax3 = Axis(fig[2, 1], 
        title="MSLP Ensemble Mean (hPa)",
        xlabel="Longitude", ylabel="Latitude")
    hm3 = heatmap!(ax3, lons, lats, plot_mslp_mean, colormap=:balance)
    Colorbar(fig[2, 2], hm3)
    
    # MSLP std
    ax4 = Axis(fig[2, 3], 
        title="MSLP Ensemble Std (hPa)",
        xlabel="Longitude", ylabel="Latitude")
    hm4 = heatmap!(ax4, lons, lats, plot_mslp_std, colormap=:viridis)
    Colorbar(fig[2, 4], hm4)
    
    # Add summary stats as text
    # Figure out which dimension corresponds to lat
    println("tas_mean shape: $(size(tas_mean)), lons: $(length(lons)), lats: $(length(lats))")
    
    # The data after time averaging should be (lon, lat) or (lat, lon)
    # Check which matches
    if size(tas_mean, 1) == length(lons) && size(tas_mean, 2) == length(lats)
        # (lon, lat) ordering
        nh_mask = lats .> 0
        tas_nh_mean = mean(tas_mean[:, nh_mask])
        mslp_nh_mean = mean(mslp_mean[:, nh_mask])
        tas_nh_std = mean(tas_std[:, nh_mask])
        mslp_nh_std = mean(mslp_std[:, nh_mask])
    else
        # (lat, lon) ordering
        nh_mask = lats .> 0
        tas_nh_mean = mean(tas_mean[nh_mask, :])
        mslp_nh_mean = mean(mslp_mean[nh_mask, :])
        tas_nh_std = mean(tas_std[nh_mask, :])
        mslp_nh_std = mean(mslp_std[nh_mask, :])
    end
    
    Label(fig[3, 1:4], 
        "NH Averages - TAS mean: $(round(tas_nh_mean, digits=2)) K, TAS spread: $(round(tas_nh_std, digits=3)) K | MSLP mean: $(round(mslp_nh_mean, digits=1)) hPa, MSLP spread: $(round(mslp_nh_std, digits=3)) hPa",
        fontsize=14)
    
    # Save to ClimaCoupler.jl directory
    output_dir = joinpath(@__DIR__, "..", "..", "..")  # ClimaCoupler.jl root
    output_file = joinpath(output_dir, "ensemble_maps_exp24_iter$(iteration).png")
    save(output_file, fig)
    println("Saved: $output_file")
    
    # Create grid layout for ensemble members (4 columns)
    n_cols = 4
    n_rows = ceil(Int, n_loaded / n_cols)
    
    # TAS ensemble members plot
    fig2 = Figure(size=(1600, 400 * n_rows))
    
    Label(fig2[0, 1:n_cols], "TAS for all $n_loaded ensemble members - Iteration $iteration", fontsize=18)
    
    for (i, tas_data) in enumerate(all_tas)
        row = div(i - 1, n_cols) + 1
        col = mod(i - 1, n_cols) + 1
        ax = Axis(fig2[row, col], title="Member $i", aspect=DataAspect())
        plot_data = needs_transpose ? permutedims(tas_data) : tas_data
        heatmap!(ax, lons, lats, plot_data, colormap=tas_cmap, colorrange=(tas_min, tas_max))
        hidedecorations!(ax)
    end
    
    # Add a single colorbar
    Colorbar(fig2[1:n_rows, n_cols+1], limits=(tas_min, tas_max), colormap=tas_cmap, label="TAS (K)")
    
    output_file2 = joinpath(output_dir, "ensemble_tas_members_exp24_iter$(iteration).png")
    save(output_file2, fig2)
    println("Saved: $output_file2")
    
    # MSLP ensemble members plot
    fig3 = Figure(size=(1600, 400 * n_rows))
    
    Label(fig3[0, 1:n_cols], "MSLP for all $n_loaded ensemble members - Iteration $iteration", fontsize=18)
    
    for (i, mslp_data) in enumerate(all_mslp)
        row = div(i - 1, n_cols) + 1
        col = mod(i - 1, n_cols) + 1
        ax = Axis(fig3[row, col], title="Member $i", aspect=DataAspect())
        plot_data = needs_transpose ? permutedims(mslp_data) : mslp_data
        heatmap!(ax, lons, lats, plot_data ./ 100, colormap=:balance, colorrange=(mslp_min, mslp_max))
        hidedecorations!(ax)
    end
    
    # Add a single colorbar
    Colorbar(fig3[1:n_rows, n_cols+1], limits=(mslp_min, mslp_max), colormap=:balance, label="MSLP (hPa)")
    
    output_file3 = joinpath(output_dir, "ensemble_mslp_members_exp24_iter$(iteration).png")
    save(output_file3, fig3)
    println("Saved: $output_file3")
    
    return fig, fig2, fig3
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    plot_ensemble_maps()
end
