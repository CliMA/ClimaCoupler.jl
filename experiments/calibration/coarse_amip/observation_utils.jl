# Utilities for processing observations/OutputVars
# Used to generate observations and compute the observation map
using ClimaAnalysis, ClimaCoupler
using Dates, LinearAlgebra, Statistics
import EnsembleKalmanProcesses as EKP

# Workaround to read ql, qi from file nicely
push!(ClimaAnalysis.Var.TIME_NAMES, "valid_time")

# Constants
const days_in_seconds = 86_400
const months = 31days_in_seconds
const years = 365days_in_seconds
const spinup_time = 3months
const start_date = DateTime(2000, 3, 1)
const first_year_start_date = DateTime(2000, 12, 1)

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/data_sources.jl"))

"""
    get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)

Return a NamedTuple of OutputVars containing all initial coarse AMIP observations.
Start date is set to `DateTime(2000, 3, 1)`. All OutputVars are resampled to the model diagnostic grid.
"""
function get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)
    diagnostic_var3d = limit_pressure_dim_to_era5_range(diagnostic_var3d)

    resample_2d(output_var) = resampled_as(output_var, diagnostic_var2d, dim_names = ["longitude", "latitude"])
    resample_3d(output_var) =
        resampled_as(output_var, diagnostic_var3d, dim_names = ["longitude", "latitude", "pressure_level"])
    resample(ov) = has_pressure(ov) ? resample_3d(ov) : resample_2d(ov)

    era5_outputvar(path) = OutputVar(path; new_start_date = start_date, shift_by = Dates.firstdayofmonth)

    rad_and_pr_obs_dict = get_obs_var_dict()

    # 2D Fields
    # TOA incoming shortwave radiation
    rsdt = resample(rad_and_pr_obs_dict["rsdt"](start_date))
    # TOA outgoing long, shortwave radiation
    rlut = resample(rad_and_pr_obs_dict["rlut"](start_date))
    rsut = resample(rad_and_pr_obs_dict["rsut"](start_date))
    # TOA clearsky outgoing long, shortwave radiation
    rsutcs = resample(rad_and_pr_obs_dict["rsutcs"](start_date))
    rlutcs = resample(rad_and_pr_obs_dict["rlutcs"](start_date))

    # TOA net radiative flux
    net_rad = rlut + rsut - rsdt
    # For some reason we need to add the start date back in
    # TODO: Make PR in climaanalysis to do this
    net_rad.attributes["start_date"] = string(start_date)

    # cloud radiative effect
    cre = rsut + rlut - rsutcs - rlutcs
    cre.attributes["start_date"] = string(start_date)

    # Precipitation
    pr = resample(rad_and_pr_obs_dict["pr"](start_date))

    # Latent heat flux
    lhf = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_averages_surface_single_level_mslhf.nc")))
    # Sensible heat flux
    shf = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_averages_surface_single_level_msshf.nc")))
    shf = lhf + shf
    shf.attributes["start_date"] = string(start_date)
    # Surface temperature
    ts = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_avg_ts.nc")))
    # 3D Fields
    # Air temperature 
    ta = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_avg_pressure_level_t.nc")))

    # relative humidity
    hur = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_avg_pressure_level_r.nc")))
    # specific humidity
    hus = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_avg_pressure_level_q.nc")))

    # # Cloud specific liquid water content
    # ql = era5_outputvar(joinpath(obs_dir, "era5_specific_cloud_liquid_water_content_1deg.nc"))
    # # Cloud specific ice water content
    # qi = era5_outputvar(joinpath(obs_dir, "era5_specific_cloud_ice_water_content_1deg.nc"))
    # foreach((ql, qi)) do var
    #     # Convert from hPa to Pa in-place so we don't create more huge OutputVars
    #     @assert var.dim_attributes[pressure_name(var)]["units"] == "hPa"
    #     var.dims[pressure_name(var)] .*= 100.0
    #     set_dim_units!(var, pressure_name(var), "Pa")
    # end
    # TODO: determine where time is spent here
    # ql = resample(reverse_dim(reverse_dim(ql, latitude_name(ql)), pressure_name(ql)))
    # qi = resample(reverse_dim(reverse_dim(qi, latitude_name(qi)), pressure_name(qi)))

    return (; rlut, rsut, rsutcs, rlutcs, pr, net_rad, cre, shf, ts, ta, hur, hus)#, ql, qi)
end

# The ERA5 pressure range is not as large as the ClimaAtmos default pressure levels,
# so we need to limit outputvars to the ERA5 pressure range (100.0 - 100_000.0 Pa) 
function limit_pressure_dim_to_era5_range(diagnostic_var3d)
    @assert has_pressure(diagnostic_var3d)
    era5_pressure_min = 100.0
    era5_pressure_max = 100_000.0
    pressure_dims = diagnostic_var3d.dims[pressure_name(diagnostic_var3d)]
    valid_pressure_levels = filter(pressure_dims) do pressure
        era5_pressure_min <= pressure <= era5_pressure_max
    end
    lowest_valid_level = minimum(valid_pressure_levels)
    highest_valid_level = maximum(valid_pressure_levels)

    return window(diagnostic_var3d, "pfull"; left = lowest_valid_level, right = highest_valid_level)
end

#####
# Processing to create EKP.ObservationSeries
#####
function to_datetime(var::OutputVar)
    start_date = DateTime(var.attributes["start_date"], dateformat"yyyy-mm-ddTHH:MM:SS")
    return [start_date + Second(t) for t in times(var)]
end
to_datetime(start_date, time) = DateTime(start_date) + Second(time)
to_datetime(time) = DateTime(start_date) + Second(time)

get_monthly_averages(simdir, var_name) = get(simdir; short_name = var_name, reduction = "average", period = "1M")

get_seasonal_averages(var) = average_time.(split_by_season_across_time(var))

function get_seasonal_averages(simdir, var_name)
    var = get(simdir; short_name = var_name, reduction = "average", period = "1M")
    get_seasonal_averages(var)
end

function get_yearly_averages(var)
    seasonal_avgs = get_seasonal_averages(var)
    nyears = fld(length(seasonal_avgs), 4)
    matrices = getproperty.(seasonal_avgs, :data)
    year_averaged_matrices = map(1:nyears) do i
        start_idx = (i - 1) * 4 + 1
        end_idx = i * 4
        group = matrices[start_idx:end_idx]

        # Compute the average matrix for this group
        averaged_matrix = sum(group) / 4
        averaged_matrix
    end
    return year_averaged_matrices
end

# TODO: compute seasonal stdev properly
function get_seasonal_stdev(output_var)
    all_seasonal_averages = get_seasonal_averages(output_var)
    all_seasonal_averages = downsample.(all_seasonal_averages, 3)
    
    # Determine dimensionality of the data
    dims = ndims(all_seasonal_averages[1]) + 1
    
    seasonal_average_matrix = cat(all_seasonal_averages...; dims)
    interannual_stdev = std(seasonal_average_matrix; dims)
    return dropdims(interannual_stdev; dims)
end

# Given an outputvar, compute the covariance for each season.
function get_seasonal_covariance(output_var)
    stdev = get_seasonal_stdev(output_var)
    return Diagonal(vec(stdev) .^ 2)
end

# Given a year, return the indices of that year within a seasonal array
# Assume each year has 4 seasons and starts at index % 4 == 1
function get_year_indices(year)
    start_index = (year * 4) - 3
    end_index = year * 4
    return start_index:end_index
end

# Take in a vector of seasonal average OutputVars and a range or single number representing the years,
# return a vector of all data within the year range
function vectorize_nyears_of_seasonal_outputvars(vec_of_vars, year_range)
    # Generate indices for all specified years
    all_year_indices = vcat(get_year_indices.(year_range)...)
    result = vcat(vec.(getproperty.(vec_of_vars[all_year_indices], :data))...)
    return result
end

# Make an EKP.Observation of a single year of seasonal averages from an OutputVar
function make_single_year_of_seasonal_observations(output_var, yr)
    seasonal_avgs = get_seasonal_averages(output_var)
    downsampled_seasonal_avg_arrays = downsample.(seasonal_avgs, 3)
    all_year_indices = vcat(get_year_indices.(yr)...)
    obs_vec = vcat(vec.(downsampled_seasonal_avg_arrays[all_year_indices])...)

    name = get(output_var.attributes, "CF_name", get(output_var.attributes, "long_name", ""))
    cov = get_seasonal_covariance(output_var)
    return EKP.Observation(obs_vec, Diagonal(repeat(cov.diag, 4)), "$(yr)_$name")
end

"""
    create_observation_vector(nt, yrs = 19)

Given a NamedTuple, produce a vector of `EKP.Observation`s, where each observation
consists of seasonally averaged fields, with the exception of globally averaged yearly radiative balance
"""
function create_observation_vector(nt, yrs = 19)
    # Starting year is 2000-12 to 2001-11
    t_start = Second(first_year_start_date - start_date).value
    rsut = window(nt.rsut, "time"; left = t_start)
    rlut = window(nt.rlut, "time"; left = t_start)
    rsutcs = window(nt.rsutcs, "time"; left = t_start)
    rlutcs = window(nt.rlutcs, "time"; left = t_start)
    cre = window(nt.cre, "time"; left = t_start)

    # Compute yearly net radiative flux separately
    net_rad = window(nt.net_rad, "time"; left = t_start) |> average_lat |> average_lon
    net_rad = get_yearly_averages(net_rad)
    net_rad_stdev = std(cat(net_rad..., dims = 3), dims = 3)
    net_rad_covariance = Diagonal(vec(net_rad_stdev) .^ 2)

    ts = window(nt.ts, "time"; left = t_start)
    pr = window(nt.pr, "time"; left = t_start)
    shf = window(nt.shf, "time"; left = t_start)

    ta = window(nt.ta, "time"; left = t_start)
    hur = window(nt.hur, "time"; left = t_start)
    hus = window(nt.hus, "time"; left = t_start)

    all_observations = map(1:yrs) do yr
        net_rad_obs = EKP.Observation(vec(net_rad[yr]), net_rad_covariance, "$(yr)_net_rad")

        rsut_obs = make_single_year_of_seasonal_observations(rsut, yr)
        rlut_obs = make_single_year_of_seasonal_observations(rlut, yr)
        rsutcs_obs = make_single_year_of_seasonal_observations(rsutcs, yr)
        rlutcs_obs = make_single_year_of_seasonal_observations(rlutcs, yr)
        cre_obs = make_single_year_of_seasonal_observations(cre, yr)
        pr_obs = make_single_year_of_seasonal_observations(pr, yr)
        # shf_obs = make_single_year_of_seasonal_observations(shf, yr)
        ts_obs = make_single_year_of_seasonal_observations(ts, yr)

        ta_obs = make_single_year_of_seasonal_observations(ta, yr)
        hur_obs = make_single_year_of_seasonal_observations(hur, yr)
        hus_obs = make_single_year_of_seasonal_observations(hus, yr)
        EKP.combine_observations([net_rad_obs, rsut_obs, rlut_obs, cre_obs, pr_obs, ts_obs, ta_obs, hur_obs, hus_obs])
    end

    return all_observations # NOT an EKP.ObservationSeries
end
# TODO: Ask  kevin to implement in ClimaAnalysis
downsample(var::ClimaAnalysis.OutputVar, n) = downsample(var.data, n)

function downsample(arr::AbstractArray, n)
    if n < 1
        error("Downsampling factor n must be at least 1.")
    end
    if ndims(arr) == 2
        return arr[1:n:end, 1:n:end]
    elseif ndims(arr) == 3
        return arr[1:n:end, 1:n:end, :]
    else
        error("Only 2D and 3D arrays are supported.")
    end
end
