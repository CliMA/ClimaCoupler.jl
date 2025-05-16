# Utilities for processing observations/OutputVars
# Used to generate observations and compute the observation map
using ClimaAnalysis, ClimaCoupler
using Dates, LinearAlgebra, Statistics
import EnsembleKalmanProcesses as EKP
import NaNStatistics

# Constants
const days_in_seconds = 86_400
const months = 31days_in_seconds
const years = 365days_in_seconds
const spinup_months = 3
const spinup_time = spinup_months * months
const start_date = DateTime(2000, 3, 1)
const first_year_start_date = DateTime(2002, 12, 1)

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/data_sources.jl"))

"""
    get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)

Return a NamedTuple of OutputVars containing all initial coarse AMIP observations.
Start date is set to `DateTime(2000, 3, 1)`. All OutputVars are resampled to the model diagnostic grid.
"""
function get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)
    # diagnostic_var3d = limit_pressure_dim_to_era5_range(diagnostic_var3d)

    resample_2d(ov) =
        resampled_as(shift_longitude(ov, -180.0, 180.0), diagnostic_var2d, dim_names = ["longitude", "latitude"])
    resample_3d(ov) = resampled_as(
        shift_longitude(ov, -180.0, 180.0),
        diagnostic_var3d,
        dim_names = ["longitude", "latitude", "pressure_level"],
    )
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
    # cloud radiative effect
    cre = rsut + rlut - rsutcs - rlutcs
    pr = resample(rad_and_pr_obs_dict["pr"](start_date))
    ta = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_avg_pressure_level_t.nc")))
    ts = resample(
        OutputVar(
            joinpath(obs_dir, "era5_monthly_avg_ts.nc");
            new_start_date = start_date,
            shift_by = Dates.firstdayofmonth,
        ),
    )

    lwp = OutputVar(joinpath(pkgdir(ClimaCoupler), "modis_lwp_iwp.nc"), "lwp", new_start_date = start_date)
    lwp = resample(lwp)
    lwp = window(lwp, "latitude"; left = -60, right = 60)
    lwp = replace(lwp, NaN => NaNStatistics.nanmean(lwp.data))

    iwp = OutputVar(joinpath(pkgdir(ClimaCoupler), "modis_lwp_iwp.nc"), "iwp", new_start_date = start_date)
    iwp = resample(iwp)
    iwp = window(iwp, "latitude"; left = -60, right = 60)
    iwp = replace(iwp, NaN => NaNStatistics.nanmean(iwp.data))

    return (; net_rad, cre, rlut, rsut, rsutcs, rlutcs, rsdt, pr, ta, ts, lwp, iwp)
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

get_monthly_averages(simdir, var_name) = get(simdir; short_name = var_name, reduction = "average", period = "1M")

function get_yearly_averages(var)
    # Get seasonal averages
    seasonal_avgs = average_time.(split_by_season_across_time(var))
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

function get_seasonal_covariance(output_var; model_error_scale = nothing, regularization = nothing)
    # Split by season
    seasons = split_by_season(output_var)

    # Calculate temporal variance of seasonal averages
    variances_per_season = map(seasons) do season
        # Average each season over its months to get one value per year
        season_across_time = split_by_season_across_time(season)
        season_across_time = filter(!isempty, season_across_time)
        seasonal_means_per_year = average_time.(season_across_time)
        dims = ["longitude", "latitude"]
        has_pressure(season) && push!(dims, "pressure")
        seasonal_means_per_year = permutedims.(seasonal_means_per_year, Ref(dims))
        # Concatenate the outputvars into a matrix
        seasonal_means_per_year_matrix = cat(getproperty.(seasonal_means_per_year, :data)..., dims = 3)
        time_mean_over_years = mean(seasonal_means_per_year_matrix, dims = 3)
        # Take variance over time of the seasonal means
        variance = var(seasonal_means_per_year_matrix, dims = 3)
        # Add model error if applicable
        if !isnothing(model_error_scale)
            @. variance += (model_error_scale * time_mean_over_years)^2
        end
        return variance
    end

    diag_cov = vcat(vec.(variances_per_season)...)

    # Add regularization if applicable
    if !isnothing(regularization)
        diag_cov .+= regularization
    end

    return Diagonal(diag_cov)
end

# Given a year, return the indices of that year within a seasonal array
# Assume each year has 4 seasons and starts at index % 4 == 1
function get_year_indices(year)
    start_index = (year * 4) - 3
    end_index = year * 4
    return start_index:end_index
end

# Todo: use dates instead of year indices
function year_of_seasonal_averages(output_var, yr)
    seasons = split_by_season_across_time(output_var)

    # Average each season over its months
    seasonal_means_per_year = average_time.(seasons)

    # Ensure dimensions are consistent for all seasons
    dims = ["longitude", "latitude"]
    has_pressure(output_var) && push!(dims, "pressure")
    seasonal_means_per_year = permutedims.(seasonal_means_per_year, Ref(dims))

    # Get data for the specific year requested
    year_ind = get_year_indices(yr)
    year_seasonal_data = seasonal_means_per_year[year_ind]
    seasons = map(x -> x.attributes["season"], year_seasonal_data)
    @assert seasons == ["DJF", "MAM", "JJA", "SON"]
    obs_vec = vcat(vec.(getproperty.(year_seasonal_data, :data))...)
    return obs_vec
end

# Make an EKP.Observation of a single year of seasonal averages from an OutputVar
function make_single_year_of_seasonal_observations(output_var, yr)
    # Split by season first to get 4 OutputVars per year
    obs_vec = year_of_seasonal_averages(output_var, yr)
    cov = get_seasonal_covariance(output_var; model_error_scale = 0.05, regularization = 1e-3)
    name = get(output_var.attributes, "CF_name", get(output_var.attributes, "long_name", ""))
    return EKP.Observation(obs_vec, cov, "$(yr)_$name")
end

"""
    create_observation_vector(nt, yrs = 17)

Given a NamedTuple, produce a vector of `EKP.Observation`s, where each observation
consists of seasonally averaged fields, with the exception of globally averaged yearly radiative balance
"""
function create_observation_vector(nt, nyears = 17)
    rsut = window(nt.rsut, "time"; left = first_year_start_date)
    rlut = window(nt.rlut, "time"; left = first_year_start_date)
    cre = window(nt.cre, "time"; left = first_year_start_date)

    ts = window(nt.ts, "time"; left = first_year_start_date)
    pr = window(nt.pr, "time"; left = first_year_start_date)
    lwp = window(nt.lwp, "time"; left = first_year_start_date)
    iwp = window(nt.iwp, "time"; left = first_year_start_date)

    # Compute yearly net radiative flux separately
    net_rad = window(nt.net_rad, "time"; left = first_year_start_date) |> average_lat |> average_lon
    net_rad = get_yearly_averages(net_rad)
    net_rad_stdev = std(cat(net_rad..., dims = 3), dims = 3)
    net_rad_covariance = Diagonal(vec(net_rad_stdev) .^ 2)

    all_observations = map(1:nyears) do yr
        @info "Creating observations for year $yr"
        net_rad_obs = EKP.Observation(vec(net_rad[yr]), net_rad_covariance, "$(yr)_net_rad")

        rsut_obs = make_single_year_of_seasonal_observations(rsut, yr)
        rlut_obs = make_single_year_of_seasonal_observations(rlut, yr)
        cre_obs = make_single_year_of_seasonal_observations(cre, yr)
        pr_obs = make_single_year_of_seasonal_observations(pr, yr)
        ts_obs = make_single_year_of_seasonal_observations(ts, yr)
        lwp_obs = make_single_year_of_seasonal_observations(lwp, yr)
        iwp_obs = make_single_year_of_seasonal_observations(iwp, yr)

        return EKP.combine_observations([net_rad_obs, cre_obs, rsut_obs, rlut_obs, pr_obs, ts_obs, lwp_obs, iwp_obs])
    end

    return all_observations # NOT an EKP.ObservationSeries
end
#= 
possible issues
- to_pressure_coordinates using wrong units (should be using Pa)
- inconsistent flattening: OutputVar dimensions ordered differently before flattening, different methods of flattening
=#
