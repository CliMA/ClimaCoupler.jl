using Statistics
import JLD2
import Dates
using Dates: Week, Month, Year, Day, Millisecond
using ClimaAnalysis
import ClimaCalibrate
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import ClimaCoupler
import EnsembleKalmanProcesses as EKP

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
    n_members = EKP.get_N_ens(ekp)
    
    # Load preprocessed observation vars to get expected dimensions
    preprocessed_vars = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/preprocessed_vars.jld2"),
    )
    
    # Get expected size from one observation variable
    # All iterations use the same observation window
    sample_date_range = first(CALIBRATE_CONFIG.sample_date_ranges)
    short_names = CALIBRATE_CONFIG.short_names
    
    # Calculate expected G dimension based on observation vector construction
    expected_dim = 0
    for short_name in short_names
        key = (short_name, sample_date_range)
        if haskey(preprocessed_vars, key)
            var = preprocessed_vars[key]
            flat_data = vec(var.data)
            valid_data = filter(!isnan, collect(skipmissing(flat_data)))
            expected_dim += length(valid_data)
        end
    end
    
    @info "Expected G dimension: $expected_dim for $n_members members"
    g_ens = fill(NaN, expected_dim, n_members)

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
    return g_ens
end

"""
    process_member_to_g_vector(diagnostics_folder_path, iteration)

Process the data of a single ensemble member and return a G vector
matching the observation vector format.

Data is normalized using the same constants as the observations for consistent EKP comparison.
"""
function process_member_to_g_vector(diagnostics_folder_path, iteration)
    short_names = CALIBRATE_CONFIG.short_names
    # All iterations use the same observation window
    sample_date_range = first(CALIBRATE_CONFIG.sample_date_ranges)
    
    # Load preprocessed observations to get target grid
    preprocessed_vars = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/preprocessed_vars.jld2"),
    )
    
    # Load normalization constants (must match what was used for observations)
    norm_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/normalization.jld2")
    if !isfile(norm_path)
        error("Normalization file not found at $norm_path. Run generate_observations.jl first.")
    end
    normalization = JLD2.load_object(norm_path)
    data_mean = normalization.mean
    data_std = normalization.std
    
    @info "Short names: $short_names, using normalization: mean=$data_mean, std=$data_std"
    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    
    all_data = Float64[]
    
    for short_name in short_names
        # Get simulation variable
        var = get_var(short_name, simdir)
        var = preprocess_var(var, sample_date_range)
        
        # Get corresponding observation variable for grid matching
        obs_key = (short_name, sample_date_range)
        if haskey(preprocessed_vars, obs_key)
            obs_var = preprocessed_vars[obs_key]
            # Resample simulation to observation grid
            obs_lons = ClimaAnalysis.longitudes(obs_var)
            obs_lats = ClimaAnalysis.latitudes(obs_var)
            var = ClimaAnalysis.resampled_as(var; lon=obs_lons, lat=obs_lats)
        end
        
        # Flatten and append, matching observation vector construction
        flat_data = vec(var.data)
        valid_data = filter(!isnan, collect(skipmissing(flat_data)))
        
        # Apply the SAME normalization as observations
        normalized_data = (valid_data .- data_mean) ./ data_std
        append!(all_data, normalized_data)
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
