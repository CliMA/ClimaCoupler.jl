using Statistics
import JLD2
# import GeoMakie
# import Makie
import Dates
using ClimaAnalysis
using DataStructures: OrderedDict
import ClimaCalibrate
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import ClimaCoupler
import ClimaCalibrate: EnsembleBuilder

include(joinpath(@__DIR__, "observation_utils.jl"))

"""
    ClimaCalibrate.observation_map(iteration)

Return G ensemble for an `iteration`.

G ensemble represents the concatenated forward model evaluations from all
ensemble members, arranged horizontally. Each individual forward model
evaluation corresponds to preprocessed, flattened simulation data from a single
ensemble member that has been matched to the corresponding observational data.
"""
function ClimaCalibrate.observation_map(iteration)
    output_dir = CALIBRATE_CONFIG.output_dir
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))

    g_ens_builder = EnsembleBuilder.GEnsembleBuilder(ekp)

    # DEBUG: Write to file for guaranteed visibility
    debug_file = joinpath(output_dir, "observation_map_debug_iter$(iteration).log")
    open(debug_file, "w") do io
        println(io, "=== observation_map for iteration $iteration ===")
        println(io, "Time: $(Dates.now())")
        println(io, "G_ens initial size: $(size(g_ens_builder.g_ens))")
        println(io, "G_ens initial non-zero: $(count(!iszero, g_ens_builder.g_ens))")
        flush(io)

        for m in 1:EKP.get_N_ens(ekp)
            member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
            simdir_path = joinpath(member_path, "amip_land/output_active")
            println(io, "\n=== Processing member $m: $simdir_path ===")
            println("Processing member $m: $simdir_path")
            flush(stdout)
            flush(io)

            try
                process_member_data!(g_ens_builder, simdir_path, m, iteration; debug_io=io)

                # Check result after each member
                col_nonzero = count(!iszero, g_ens_builder.g_ens[:, m])
                col_mean = Statistics.mean(g_ens_builder.g_ens[:, m])
                println(io, "After member $m: non-zero=$(col_nonzero), mean=$(col_mean)")
                println("After member $m: non-zero=$(col_nonzero), mean=$(col_mean)")
                flush(stdout)
                flush(io)
            catch e
                println(io, "ERROR for member $m: $e")
                println(io, catch_backtrace())
                @error "Ensemble member $m failed" exception = (e, catch_backtrace())
                EnsembleBuilder.fill_g_ens_col!(g_ens_builder, m, NaN)
                flush(io)
            end
        end

        # Final stats
        println(io, "\n=== Final G_ens stats ===")
        println(io, "non-zero: $(count(!iszero, g_ens_builder.g_ens))")
        println(io, "mean: $(Statistics.mean(g_ens_builder.g_ens))")
        println(io, "range: $(extrema(g_ens_builder.g_ens))")
    end

    if count(isnan, g_ens_builder.g_ens) > 0.9 * length(g_ens_builder.g_ens)
        error("Too many NaNs")
    end
    return g_ens_builder.g_ens
end

"""
    process_member_data!(diagnostics_folder_path, short_names, current_minibatch)

Process the data of a single ensemble member and return a single column of the
G ensemble matrix.
"""
function process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx, iteration; debug_io=nothing)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges[iteration + 1]

    function debug_print(msg)
        println(msg)
        flush(stdout)
        if !isnothing(debug_io)
            println(debug_io, msg)
            flush(debug_io)
        end
    end

    debug_print("Short names: $short_names")

    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in short_names
        var = get_var(short_name, simdir)
        var = preprocess_var(var, sample_date_ranges)

        # DEBUG: Print preprocessed var info
        debug_print("=== DEBUG: Preprocessed var for $short_name ===")
        debug_print("  var.dims keys: $(keys(var.dims))")
        for (k, v) in var.dims
            if k == "time"
                debug_print("    time: $v (type: $(typeof(v)), eltype: $(eltype(v)))")
            else
                debug_print("    $k: length=$(length(v)), range=$(extrema(v))")
            end
        end
        var_short_name_attr = get(var.attributes, "short_name", "MISSING")
        var_start_date_attr = get(var.attributes, "start_date", "MISSING")
        debug_print("  var.attributes short_name: $var_short_name_attr")
        debug_print("  var.attributes start_date: $var_start_date_attr")
        debug_print("  var.data size: $(size(var.data))")
        debug_print("  var.data range: $(extrema(var.data))")
        debug_print("  var.data mean: $(Statistics.mean(var.data))")

        # DEBUG: Print metadata info from g_ens_builder
        metadata_infos = get(g_ens_builder.metadata_by_short_name, short_name, [])
        debug_print("  Found $(length(metadata_infos)) metadata entries for $short_name")
        for (i, md_info) in enumerate(metadata_infos)
            md = md_info.metadata
            debug_print("  Metadata $i (range=$(md_info.range)):")
            debug_print("    dims keys: $(keys(md.dims))")
            for (k, v) in md.dims
                if k == "time"
                    debug_print("      time: $v (type: $(typeof(v)), eltype: $(eltype(v)))")
                else
                    debug_print("      $k: length=$(length(v)), range=$(extrema(v))")
                end
            end
            md_start_date = get(md.attributes, "start_date", "MISSING")
            debug_print("    start_date: $md_start_date")
        end

        result = EnsembleBuilder.fill_g_ens_col!(g_ens_builder, col_idx, var; verbose = true)
        debug_print("  fill_g_ens_col! returned: $result")

        # Check what's in the g_ens column now
        col_data = g_ens_builder.g_ens[:, col_idx]
        debug_print("  After fill: col non-zero=$(count(!iszero, col_data)), mean=$(Statistics.mean(col_data))")
    end

    return nothing
end

function largest_period(sample_date_range)
    span = maximum(sample_date_range) - minimum(sample_date_range)
    span = Millisecond(span)
    span.value == 0 && return Month(1)
    day_in_ms = 8.64e7
    period =
        span.value >= day_in_ms * 365 ? Year(1) :
        span.value >= day_in_ms * 30 ? Month(1) :
        span.value >= day_in_ms * 7 ? Week(1) : Day(1)
    return period
end

# 3D variables that need z-to-pressure conversion
const PRESSURE_LEVEL_VARS = ["ta", "ua", "va", "hus", "hur"]

# Standard ERA5 pressure levels in hPa (must match generate_observations.jl)
# Note: 1000 and 925 hPa excluded because model z-windowing (left=80) removes near-surface levels
const ERA5_PRESSURE_LEVELS = [850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 250.0, 200.0]

function get_var(short_name, simdir)
    if short_name == "tas - ta"
        tas = get(simdir; short_name = "tas", period = "1M")
        ta = get(simdir; short_name = "ta", period = "1M")
        # TODO: figure out why this doesn't work
        # pfull = get(simdir; short_name = "pfull", period = "1M")
        # ta_900 = ClimaAnalysis.Atmos.to_pressure_coordinates(ta, pfull; target_pressure=[900])
        ta_900hpa = slice(ta; z = 1000)
        var = tas - ta_900hpa
    else
        # Specify period to avoid ambiguity when multiple periods exist
        var = get(simdir; short_name, period = "1M")
    end

    # Convert 3D variables from z-levels to pressure coordinates
    if short_name in PRESSURE_LEVEL_VARS && ClimaAnalysis.has_altitude(var)
        pfull = get(simdir; short_name = "pfull", period = "1M")
        # Window to avoid pressure inversions at low elevation
        pfull_windowed = ClimaAnalysis.window(pfull, "z", left = 80)
        var_windowed = ClimaAnalysis.window(var, "z", left = 80)
        var = ClimaAnalysis.Atmos.to_pressure_coordinates(var_windowed, pfull_windowed)
        var = ClimaAnalysis.Var.convert_dim_units(
            var,
            "pfull",
            "hPa";
            conversion_function = x -> 0.01 * x,
        )
        # Note: pressure level resampling is done in preprocess_var after time handling
    end

    var.attributes["short_name"] = short_name
    return var
end

"""
    preprocess_var(var::ClimaAnalysis.OutputVar, reference_date)

Preprocess `var` before flattening for G ensemble matrix.

For "pr", weekly sums are computed. For "tas" and "mslp", weekly means are
computed from daily means. The daily means are computing starting from
`reference_date`.

This function assumes that the data is monthly.
"""
function preprocess_var(var, sample_date_range)
    period = largest_period(sample_date_range)
    var = ClimaAnalysis.Var._shift_by(var, date -> date - period)
    var = set_units(var, var_units[short_name(var)])
    # TODO: Match dates instead of just windowing
    var = window(var, "time"; left = sample_date_range[1], right = sample_date_range[2])

    # Note: Do NOT average over time here - the time dimension must be preserved
    # to match the observation metadata dimensions. The EnsembleBuilder expects
    # simulation and observation data to have matching dimension names.

    # Resample 3D pressure-level variables to ERA5 standard levels
    # Must handle time dimension specially - slice it out, resample 3D, then add back
    if ClimaAnalysis.has_pressure(var)
        if ClimaAnalysis.has_time(var)
            # Get the time value(s) before slicing
            time_vals = ClimaAnalysis.times(var)
            # For single time point, slice to remove time dim, resample, then wrap back
            if length(time_vals) == 1
                time_val = time_vals[1]
                # Slice to get 3D var (removes time dimension)
                var_3d = ClimaAnalysis.slice(var; time = time_val)
                # Resample pressure on 3D data
                lon = ClimaAnalysis.longitudes(var_3d)
                lat = ClimaAnalysis.latitudes(var_3d)
                var_resampled = ClimaAnalysis.resampled_as(var_3d; lon, lat, pressure_level = ERA5_PRESSURE_LEVELS)

                # Rename dimensions to match ERA5 observation metadata
                # resampled_as keeps original dim names, but we need to match metadata
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

                # Add time dimension back by creating new OutputVar with time axis
                # Use merge to preserve OrderedDict types that OutputVar expects
                time_dim = OrderedDict("time" => [time_val])
                new_dims = merge(renamed_dims, time_dim)
                # Also rename dim_attributes keys to match
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
                time_dim_attribs = OrderedDict("time" => Dict{String, Any}())
                new_dim_attribs = merge(renamed_dim_attribs, time_dim_attribs)
                # Reshape data to add time dimension (last axis)
                new_data = reshape(var_resampled.data, size(var_resampled.data)..., 1)
                # Use 2-arg constructor which handles type inference better
                var = ClimaAnalysis.OutputVar(new_dims, new_data)
                # Copy over attributes
                merge!(var.attributes, var_resampled.attributes)
                merge!(var.dim_attributes, new_dim_attribs)
            else
                # Multiple time points - process each separately and concatenate
                error("Multiple time points not yet supported in preprocess_var")
            end
        else
            # No time dimension - just resample directly
            lon = ClimaAnalysis.longitudes(var)
            lat = ClimaAnalysis.latitudes(var)
            var_resampled = ClimaAnalysis.resampled_as(var; lon, lat, pressure_level = ERA5_PRESSURE_LEVELS)

            # Rename dimensions to match ERA5 observation metadata
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
            var = ClimaAnalysis.OutputVar(renamed_dims, var_resampled.data)
            merge!(var.attributes, var_resampled.attributes)
            merge!(var.dim_attributes, renamed_dim_attribs)
        end
    end
    return var
end

"""
    ClimaCalibrate.analyze_iteration(ekp,
                                     g_ensemble,
                                     prior,
                                     output_dir,
                                     iteration)

Analyze an iteration by plotting the bias plots, constrained parameters over
iterations, and errors over iterations and time.
"""
function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
    plot_constrained_params_and_errors(plot_output_path, ekp, prior)

    simdir = SimDir(ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1))
    plot_bias(ekp, simdir, iteration; output_dir = plot_output_path)

    @info "Ensemble spread: $(scalar_spread(ekp))"
end

"""
    plot_constrained_params_and_errors(output_dir, ekp, prior)

Plot the constrained parameters and errors from `ekp` and `prior` and save
them to `output_dir`.
"""
function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 1], ekp, error_metric = "loss")
    EKP.Visualize.plot_error_over_time(fig[1, dim_size + 2], ekp, error_metric = "loss")
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
end

bias_plot_extrema = Dict(
    "tas" => (-6, 6),
    "tas - ta" => (-6, 6),
    "hfls" => (-50, 50),
    "hfss" => (-25, 25),
    "rsus" => (-50, 50),
    "rlus" => (-50, 50),
    "mslp" => (-1000, 1000),
    "pr" => (-1e-4, 1e-4),
    "ta" => (-10, 10),
    "ua" => (-10, 10),
    "va" => (-5, 5),
)

"""
    plot_bias(ekp, simdir, iteration; output_dir)

Plot bias maps comparing simulation output to observations for all variables in
`CALIBRATE_CONFIG.short_names`.

Uses observations from the EKP object via `reconstruct_vars` and compares
them against simulation variables for each date in the sample date ranges.
"""
function plot_bias(ekp, simdir, iteration; output_dir = simdir.simulation_path)
    # Get observations for this iteration
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

    # Get simulation variables for all short_names
    # Note: Don't use preprocess_var here as it removes the time dimension
    # which is needed for slicing by date
    sim_vars = []
    for short_name in CALIBRATE_CONFIG.short_names
        var = get_var(short_name, simdir)
        # Apply preprocessing steps except time averaging
        period = largest_period(sample_date_ranges)
        var = ClimaAnalysis.Var._shift_by(var, date -> date - period)
        var = set_units(var, var_units[short_name])
        var = window(var, "time"; left = sample_date_ranges[1], right = sample_date_ranges[2])
        # Resample 3D pressure-level variables to ERA5 standard levels
        # Must handle time dimension - slice it out, resample 3D, then add back
        if ClimaAnalysis.has_pressure(var)
            if ClimaAnalysis.has_time(var)
                time_vals = ClimaAnalysis.times(var)
                if length(time_vals) == 1
                    time_val = time_vals[1]
                    var_3d = ClimaAnalysis.slice(var; time = time_val)
                    lon = ClimaAnalysis.longitudes(var_3d)
                    lat = ClimaAnalysis.latitudes(var_3d)
                    var_resampled = ClimaAnalysis.resampled_as(var_3d; lon, lat, pressure_level = ERA5_PRESSURE_LEVELS)

                    # Rename dimensions to match ERA5 observation metadata
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

                    # Add time dimension back using merge to preserve types
                    time_dim = OrderedDict("time" => [time_val])
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
                    time_dim_attribs = OrderedDict("time" => Dict{String, Any}())
                    new_dim_attribs = merge(renamed_dim_attribs, time_dim_attribs)
                    new_data = reshape(var_resampled.data, size(var_resampled.data)..., 1)
                    var = ClimaAnalysis.OutputVar(new_dims, new_data)
                    merge!(var.attributes, var_resampled.attributes)
                    merge!(var.dim_attributes, new_dim_attribs)
                else
                    error("Multiple time points not yet supported in plot_bias preprocessing")
                end
            else
                lon = ClimaAnalysis.longitudes(var)
                lat = ClimaAnalysis.latitudes(var)
                var_resampled = ClimaAnalysis.resampled_as(var; lon, lat, pressure_level = ERA5_PRESSURE_LEVELS)

                # Rename dimensions to match ERA5 observation metadata
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
                var = ClimaAnalysis.OutputVar(renamed_dims, var_resampled.data)
                merge!(var.attributes, var_resampled.attributes)
                merge!(var.dim_attributes, renamed_dim_attribs)
            end
        end
        push!(sim_vars, var)
    end

    # Match sim_vars with era5_vars by short_name
    var_pairs = []
    for sim_var in sim_vars
        sim_short_name = ClimaAnalysis.short_name(sim_var)
        era5_idx = findfirst(v -> ClimaAnalysis.short_name(v) == sim_short_name, era5_vars)
        if !isnothing(era5_idx)
            push!(var_pairs, (sim_var, era5_vars[era5_idx]))
        else
            @warn "No ERA5 data found for $(sim_short_name)"
        end
    end

    if isempty(var_pairs)
        @warn "No matching variable pairs found for bias plotting"
        return nothing
    end

    # Get sample dates for this iteration
    sample_dates = unique(CALIBRATE_CONFIG.sample_date_ranges[iteration + 1])

    # Use 500 hPa as representative level for bias plots
    PLOT_PRESSURE_LEVEL = 500.0

    # Create figure
    fig = GeoMakie.Figure(size = (1500, 500 * length(var_pairs)))
    for (j, date) in enumerate(sample_dates)
        for (i, (sim_var, era5_var)) in enumerate(var_pairs)
            sim_var_t = slice(sim_var, time = date)
            era5_var_t = slice(era5_var, time = date)

            # For 3D pressure-level data, slice at representative level to get 2D
            if ClimaAnalysis.has_pressure(sim_var_t)
                sim_var_t = slice(sim_var_t; pressure_level = PLOT_PRESSURE_LEVEL)
                era5_var_t = slice(era5_var_t; pressure_level = PLOT_PRESSURE_LEVEL)
            end

            # Calculate biases
            global_bias = ClimaAnalysis.global_bias(sim_var_t, era5_var_t)
            global_mean = weighted_average_lonlat(sim_var_t).data[1]
            relative_global_bias = global_bias / global_mean

            land_bias = ClimaAnalysis.global_bias(
                sim_var_t,
                era5_var_t;
                mask = ClimaAnalysis.apply_oceanmask,
            )
            land_mean =
                weighted_average_lonlat(ClimaAnalysis.apply_oceanmask(sim_var_t)).data[1]
            relative_land_bias = land_bias / land_mean

            ocean_bias = ClimaAnalysis.global_bias(
                sim_var_t,
                era5_var_t;
                mask = ClimaAnalysis.apply_landmask,
            )
            ocean_mean =
                weighted_average_lonlat(ClimaAnalysis.apply_landmask(sim_var_t)).data[1]
            relative_ocean_bias = ocean_bias / ocean_mean

            @info short_name(sim_var_t) relative_global_bias global_bias global_mean
            @info short_name(sim_var_t) relative_land_bias land_bias land_mean
            @info short_name(sim_var_t) relative_ocean_bias ocean_bias ocean_mean

            cmap_extrema =
                get(bias_plot_extrema, short_name(sim_var_t), extrema(sim_var_t.data))

            # Plot bias
            ax = ClimaAnalysis.Visualize.plot_bias_on_globe!(
                fig[i, j],
                sim_var_t,
                era5_var_t;
                cmap_extrema,
            )
        end
    end

    GeoMakie.save(joinpath(output_dir, "bias_sample_dates.png"), fig)
    return nothing
end


function scalar_spread(ekp)
    g_mean_final = EKP.get_g_mean_final(ekp)
    g_final = EKP.get_g_final(ekp)
    sq_dists = [sum((col .- g_mean_final) .^ 2) for col in eachcol(g_final)]
    return mean(sq_dists)
end
