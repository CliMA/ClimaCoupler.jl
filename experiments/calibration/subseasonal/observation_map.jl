using Statistics
import JLD2
import Dates
using Dates: Week, Month, Year, Day, Millisecond
using ClimaAnalysis
import ClimaCalibrate
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import ClimaCoupler
import EnsembleKalmanProcesses as EKP

# Override JLD2's default_iotype to use IOStream instead of MmapIO
# This avoids Bus errors from memory-mapped files on Lustre filesystem
JLD2.default_iotype() = IOStream

include(joinpath(@__DIR__, "observation_utils.jl"))

"""
    retry_on_error(f; max_retries=5, delay=2.0, backoff=2.0)

Retry a function `f` on any error, with exponential backoff.
Useful for handling transient filesystem errors on Lustre (NetCDF/HDF5 race conditions).
"""
function retry_on_error(f; max_retries=5, delay=2.0, backoff=2.0)
    last_error = nothing
    for attempt in 1:max_retries
        try
            return f()
        catch e
            last_error = e
            if attempt < max_retries
                wait_time = delay * backoff^(attempt - 1)
                @warn "Attempt $attempt failed, retrying in $(wait_time)s..." exception=(e, catch_backtrace())
                sleep(wait_time)
            end
        end
    end
    @error "All $max_retries attempts failed"
    throw(last_error)  # Use throw() not rethrow() outside catch block
end

"""
    ClimaCalibrate.observation_map(iteration)

Return G ensemble for an `iteration`.

G ensemble represents the concatenated forward model evaluations from all
ensemble members, arranged horizontally. Each individual forward model
evaluation corresponds to preprocessed, flattened simulation data from a single
ensemble member that has been matched to the corresponding observational data.

**MODES** (controlled by use_nh_averages saved in land_mask.jld2):
1. Full gridpoint mode: returns 113k+ values per member
2. NH averages mode: returns 3 values per member (one NH average per variable)

Land masking is applied per-variable:
- Variables in land_only_vars (tas, pr): only land points included
- Other variables (mslp): all points included (evaluated everywhere)
"""
function ClimaCalibrate.observation_map(iteration)
    output_dir = CALIBRATE_CONFIG.output_dir
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
    n_members = EKP.get_N_ens(ekp)
    
    # Load preprocessed observation vars to get expected dimensions
    preprocessed_vars = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/preprocessed_vars.jld2"),
    )
    
    # Load land mask, land_only_vars, and mode info
    land_mask_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/land_mask.jld2")
    land_mask = nothing
    land_only_vars = Set{String}()
    use_nh_averages = false
    if isfile(land_mask_path)
        land_mask_data = JLD2.load_object(land_mask_path)
        land_mask = land_mask_data.mask
        # Load land_only_vars if present, otherwise default to ["tas", "pr"]
        land_only_vars = Set(get(land_mask_data, :land_only_vars, ["pr"]))
        # Load NH averages mode flag
        use_nh_averages = get(land_mask_data, :use_nh_averages, false)
        @info "Loaded land mask: $(count(land_mask)) land points, land-only vars: $land_only_vars"
        @info "Observation mode: $(use_nh_averages ? "NH averages" : "full gridpoints")"
    end
    
    # Get expected size from one observation variable
    # All iterations use the same observation window
    sample_date_range = first(CALIBRATE_CONFIG.sample_date_ranges)
    short_names = CALIBRATE_CONFIG.short_names
    
    # Calculate expected G dimension based on mode
    if use_nh_averages
        # NH averages mode: one value per variable
        expected_dim = length(short_names)
        @info "NH averages mode: $expected_dim values (one per variable)"
    else
        # Full gridpoint mode: calculate based on observation vector construction
        expected_dim = 0
        for short_name in short_names
            key = (short_name, sample_date_range)
            if haskey(preprocessed_vars, key)
                var = preprocessed_vars[key]
                lats = ClimaAnalysis.latitudes(var)
                nlon = length(ClimaAnalysis.longitudes(var))
                nlat = length(lats)

                # Build combined mask: land + latitude >= -50
                combined_mask = if !isnothing(land_mask) && short_name in land_only_vars
                    copy(land_mask)
                else
                    trues(nlon, nlat)
                end
                for j in 1:nlat
                    if lats[j] < -90.0
                        combined_mask[:, j] .= false
                    end
                end

                flat_data = vec(var.data)[findall(vec(combined_mask))]
                valid_data = filter(!isnan, collect(skipmissing(flat_data)))
                expected_dim += length(valid_data)
                
                mask_status = (!isnothing(land_mask) && short_name in land_only_vars) ? "land only" : "all points"
                @info "  $short_name: $(length(valid_data)) points ($mask_status, lat >= -50)"
            end
        end
    end
    
    @info "Expected G dimension: $expected_dim for $n_members members"
    g_ens = fill(NaN, expected_dim, n_members)

    # Small delay to ensure all files are flushed to Lustre before reading
    sleep(2.0)
    
    for m in 1:n_members
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        # Find the job output folder (named after job_id, e.g., "wxquest_diagedmf")
        # It's the first directory in member_path that isn't a file
        subdirs = filter(isdir, [joinpath(member_path, f) for f in readdir(member_path)])
        if isempty(subdirs)
            error("No output directory found in $member_path")
        end
        job_dir = first(subdirs)
        # SimDir expects the clima_atmos subdirectory containing NetCDF files
        simdir_path = joinpath(job_dir, "output_active", "clima_atmos")
        @info "Processing member $m: $simdir_path"
        try
            g_col = process_member_to_g_vector(simdir_path, iteration)
            g_ens[:, m] = g_col
        catch e
            @error "Ensemble member $m failed" exception = (e, catch_backtrace())
            g_ens[:, m] .= NaN
        end
    end
    
    if count(isnan, g_ens) > 0.9 * length(g_ens)
        error("Too many NaNs in G ensemble")
    end
    
    # ==========================================================================
    # DIAGNOSTIC: Compare G ensemble mean vs y (observations)
    # This helps verify that model and observations are properly aligned
    # ==========================================================================
    # obs_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/obs_vec.jld2")
    # if isfile(obs_path)
    #     obs_vec = JLD2.load_object(obs_path)
    #     if !isempty(obs_vec)
    #         y = first(first(obs_vec).samples)  # Get the observation vector (samples is a Vector of Vectors)
    #         g_mean = vec(Statistics.mean(g_ens, dims=2))
            
    #         if length(y) == length(g_mean)
    #             @info "=== DIAGNOSTIC: G vs y comparison (normalized units) ==="
                
    #             if use_nh_averages
    #                 # NH averages mode: simple per-variable comparison (one value each)
    #                 for (i, short_name) in enumerate(short_names)
    #                     bias = g_mean[i] - y[i]
    #                     @info "  $short_name: G=$(round(g_mean[i], sigdigits=4)), y=$(round(y[i], sigdigits=4)), bias=$(round(bias, sigdigits=3))"
    #                 end
    #                 overall_rmse = sqrt(Statistics.mean((g_mean .- y).^2))
    #                 @info "  Overall RMSE: $(round(overall_rmse, sigdigits=3))"
    #             else
    #                 # Full gridpoint mode: compute per-variable statistics
    #                 tas_end = 24472
    #                 mslp_end = tas_end + 64800
    #                 pr_end = mslp_end + 24472
                    
    #                 if length(y) >= pr_end
    #                     # Per-variable bias (G - y)
    #                     tas_bias = Statistics.mean(g_mean[1:tas_end] .- y[1:tas_end])
    #                     mslp_bias = Statistics.mean(g_mean[tas_end+1:mslp_end] .- y[tas_end+1:mslp_end])
    #                     pr_bias = Statistics.mean(g_mean[mslp_end+1:pr_end] .- y[mslp_end+1:pr_end])
                        
    #                     # Per-variable RMSE
    #                     tas_rmse = sqrt(Statistics.mean((g_mean[1:tas_end] .- y[1:tas_end]).^2))
    #                     mslp_rmse = sqrt(Statistics.mean((g_mean[tas_end+1:mslp_end] .- y[tas_end+1:mslp_end]).^2))
    #                     pr_rmse = sqrt(Statistics.mean((g_mean[mslp_end+1:pr_end] .- y[mslp_end+1:pr_end]).^2))
                        
    #                     @info "  tas:  bias=$(round(tas_bias, sigdigits=3)), RMSE=$(round(tas_rmse, sigdigits=3))"
    #                     @info "  mslp: bias=$(round(mslp_bias, sigdigits=3)), RMSE=$(round(mslp_rmse, sigdigits=3))"
    #                     @info "  pr:   bias=$(round(pr_bias, sigdigits=3)), RMSE=$(round(pr_rmse, sigdigits=3))"
    #                 end
    #                 @info "  Overall RMSE: $(round(sqrt(Statistics.mean((g_mean .- y).^2)), sigdigits=3))"
    #             end
    #         else
    #             @warn "Dimension mismatch: y=$(length(y)), G=$(length(g_mean))"
    #         end
    #     end
    # end
    
    return g_ens
end

"""
    process_member_to_g_vector(diagnostics_folder_path, iteration)

Process the data of a single ensemble member and return a G vector
matching the observation vector format.

**MODES** (controlled by use_nh_averages saved in land_mask.jld2):
1. Full gridpoint mode: returns 113k+ values
2. NH averages mode: returns 3 values (one NH average per variable)

**PER-VARIABLE NORMALIZATION**: Each variable is normalized using the SAME per-variable
constants computed from observations. This ensures model output is compared consistently
with observations, with each variable contributing equally to the objective.

Land masking is applied per-variable to match observation vector construction:
- Variables in land_only_vars (tas, pr): only land points included
- Other variables (mslp): all points included
"""
function process_member_to_g_vector(diagnostics_folder_path, iteration)
    short_names = CALIBRATE_CONFIG.short_names
    # All iterations use the same observation window
    sample_date_range = first(CALIBRATE_CONFIG.sample_date_ranges)
    
    # Load preprocessed observations to get target grid
    preprocessed_vars = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/preprocessed_vars.jld2"),
    )
    
    # Load per-variable normalization constants (must match what was used for observations)
    norm_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/normalization.jld2")
    if !isfile(norm_path)
        error("Normalization file not found at $norm_path. Run generate_observations.jl first.")
    end
    var_normalization = JLD2.load_object(norm_path)
    
    # Load land mask, land_only_vars, and mode info
    land_mask_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/land_mask.jld2")
    land_mask = nothing
    land_only_vars = Set{String}()
    use_nh_averages = false
    if isfile(land_mask_path)
        land_mask_data = JLD2.load_object(land_mask_path)
        land_mask = land_mask_data.mask
        # Load land_only_vars if present, otherwise default to ["tas", "pr"]
        land_only_vars = Set(get(land_mask_data, :land_only_vars, ["pr"]))
        # Load NH averages mode flag
        use_nh_averages = get(land_mask_data, :use_nh_averages, false)
    end
    
    @info "Short names: $short_names, per-variable normalization, land_only_vars=$land_only_vars"

    # Wrap SimDir in retry logic to handle transient Lustre/NetCDF errors
    simdir = retry_on_error() do
        ClimaAnalysis.SimDir(diagnostics_folder_path)
    end
    
    all_data = Float64[]
    
    for short_name in short_names
        # Get simulation variable (with retry for NetCDF access)
        var = retry_on_error() do
            get_var(short_name, simdir)
        end
        var = preprocess_var(var, sample_date_range)
        
        # Get corresponding observation variable for grid matching
        obs_key = (short_name, sample_date_range)
        if haskey(preprocessed_vars, obs_key)
            obs_var = preprocessed_vars[obs_key]
            # Resample simulation to observation grid
            obs_lons = ClimaAnalysis.longitudes(obs_var)
            obs_lats = ClimaAnalysis.latitudes(obs_var)
            
            # DIAGNOSTIC: Print grid info on first variable to verify axis ordering
            if short_name == first(short_names)
                model_lons = ClimaAnalysis.longitudes(var)
                model_lats = ClimaAnalysis.latitudes(var)
                @info "Grid check: model=$(size(var.data)), obs=$(size(obs_var.data))"
                @info "  Model lons: [$(model_lons[1]), ..., $(model_lons[end])], n=$(length(model_lons))"
                @info "  Model lats: [$(model_lats[1]), ..., $(model_lats[end])], n=$(length(model_lats))"
                @info "  Obs lons: [$(obs_lons[1]), ..., $(obs_lons[end])], n=$(length(obs_lons))"
                @info "  Obs lats: [$(obs_lats[1]), ..., $(obs_lats[end])], n=$(length(obs_lats))"
            end
            
            var = ClimaAnalysis.resampled_as(var; lon=obs_lons, lat=obs_lats)
        end
        
        # Apply precipitation unit conversion for model data
        # Model stores pr as NEGATIVE kg/m²/s, convert to positive mm/day
        if short_name == "pr"
            raw_mean = Statistics.mean(filter(!isnan, vec(var.data)))
            converted_data = var.data .* MODEL_PR_CONVERSION
            conv_mean = Statistics.mean(filter(!isnan, vec(converted_data)))
            @info "  pr unit conversion: $(round(raw_mean, sigdigits=3)) kg/m²/s → $(round(conv_mean, sigdigits=3)) mm/day"
            var = ClimaAnalysis.OutputVar(var.attributes, var.dims, var.dim_attributes, converted_data)
        end
        
        # Mode-specific processing
        if use_nh_averages
            # NH AVERAGES MODE: compute single NH average for this variable
            use_land = short_name in land_only_vars
            nh_avg = compute_nh_average(var, land_mask, use_land)
            
            # Normalize using the same constants as observations
            norm = var_normalization[short_name]
            normalized_val = (nh_avg - norm.mean) / norm.std
            push!(all_data, normalized_val)
            
            mask_status = use_land ? "NH land" : "NH all"
            @info "  $short_name: NH_avg=$(round(nh_avg, sigdigits=6)) → normalized=$(round(normalized_val, sigdigits=4)) ($mask_status)"
        else
            # FULL GRIDPOINT MODE: flatten and append
            lats = ClimaAnalysis.latitudes(var)
            nlon = length(ClimaAnalysis.longitudes(var))
            nlat = length(lats)

            # Build combined mask: land + latitude >= -50 (MUST match observation vector construction)
            combined_mask = if !isnothing(land_mask) && short_name in land_only_vars
                copy(land_mask)
            else
                trues(nlon, nlat)
            end
            for j in 1:nlat
                if lats[j] < -90.0
                    combined_mask[:, j] .= false
                end
            end

            flat_data = vec(var.data)[findall(vec(combined_mask))]
            valid_data = filter(!isnan, collect(skipmissing(flat_data)))
            
            # Apply PER-VARIABLE normalization (MUST match observation vector construction)
            norm = var_normalization[short_name]
            normalized_data = (valid_data .- norm.mean) ./ norm.std
            append!(all_data, normalized_data)
            
            mask_status = (!isnothing(land_mask) && short_name in land_only_vars) ? "land only" : "all points"
            @info "  $short_name: $(length(valid_data)) points ($mask_status, lat >= -50)"
            
            # Diagnostic: print statistics for debugging alignment issues
            @info "    G stats: mean=$(round(Statistics.mean(normalized_data), sigdigits=4)), std=$(round(Statistics.std(normalized_data), sigdigits=4)), range=[$(round(minimum(normalized_data), sigdigits=3)), $(round(maximum(normalized_data), sigdigits=3))]"
        end
    end

    return all_data
end

function largest_period(sample_date_range)
    span = maximum(sample_date_range) - minimum(sample_date_range)
    span = Millisecond(span)
    # For weekly data (6 days span = 7 day period including endpoints)
    span.value == 0 && return Week(1)
    day_in_ms = 8.64e7
    period =
        span.value >= day_in_ms * 365 ? Year(1) :
        span.value >= day_in_ms * 28 ? Month(1) :
        span.value >= day_in_ms * 6 ? Week(1) : Day(1)
    return period
end

function get_var(short_name, simdir)
    if short_name == "tas - ta"
        tas = get(simdir; short_name = "tas")
        ta = get(simdir; short_name = "ta")
        ta_900hpa = slice(ta; z = 1000)
        var = tas - ta_900hpa
    else
        var = get(simdir; short_name)
    end

    var.attributes["short_name"] = short_name
    return var
end

"""
    preprocess_var(var::ClimaAnalysis.OutputVar, sample_date_range)

Preprocess `var` before flattening for G ensemble matrix.
For weekly data, compute the temporal average over the sample window.
"""
function preprocess_var(var, sample_date_range)
    # If var has time dimension, handle it
    if "time" in keys(var.dims)
        n_times = length(ClimaAnalysis.times(var))
        if n_times > 1
            # select the time window for model output
            var = window(var, "time"; left = sample_date_range[1], right = sample_date_range[2])
            # Multiple time steps - average over time
            var = ClimaAnalysis.average_time(var)
        else
            # Single time step - just slice to remove the dimension
            var = ClimaAnalysis.slice(var; time = first(ClimaAnalysis.times(var)))
        end
    end
    
    var = set_units(var, var_units[short_name(var)])
    return var
end

"""
    ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)

Placeholder for iteration analysis (bias plots, etc.)
"""
function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    @info "Iteration $iteration analysis (placeholder)"
    # TODO: Add bias plotting when needed
    return nothing
end
