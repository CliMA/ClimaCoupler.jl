# Utilities for processing observations/OutputVars
# Used to generate observations and compute the observation map
using ClimaAnalysis, ClimaCoupler
import ClimaCalibrate
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

set_short_name!(var, short_name) = (var.attributes["short_name"] = short_name)

"""
    get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)

Return a NamedTuple of OutputVars containing all initial coarse AMIP observations.
Start date is set to `DateTime(2000, 3, 1)`. All OutputVars are resampled to the model diagnostic grid.
"""
function get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)
    if !isnothing(diagnostic_var3d)
        diagnostic_var3d = limit_pressure_dim_to_era5_range(diagnostic_var3d)
    end

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
    set_short_name!(rlut, "rlut")
    rsut = resample(rad_and_pr_obs_dict["rsut"](start_date))
    set_short_name!(rsut, "rsut")
    # TOA clearsky outgoing long, shortwave radiation
    rsutcs = resample(rad_and_pr_obs_dict["rsutcs"](start_date))
    set_short_name!(rsutcs, "rsutcs")
    rlutcs = resample(rad_and_pr_obs_dict["rlutcs"](start_date))
    set_short_name!(rlutcs, "rlutcs")

    # TOA net radiative flux
    net_rad = rlut + rsut - rsdt |> average_lonlat
    set_short_name!(net_rad, "net_rad")
    # cloud radiative effect
    sw_cre = rsut - rsutcs
    set_short_name!(sw_cre, "sw_cre")
    lw_cre = rlut - rlutcs
    set_short_name!(lw_cre, "lw_cre")
    pr = resample(rad_and_pr_obs_dict["pr"](start_date))
    set_short_name!(pr, "pr")
    # ta = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_avg_pressure_level_t.nc")))
    ts = resample(
        OutputVar(
            joinpath(obs_dir, "era5_monthly_avg_ts.nc");
            new_start_date = start_date,
            shift_by = Dates.firstdayofmonth,
        ),
    )
    set_short_name!(ts, "ts")

    lwp, iwp = map(["lwp", "iwp"]) do short_name
        out_var = OutputVar(joinpath(pkgdir(ClimaCoupler), "modis_lwp_iwp.nc"), short_name, new_start_date = start_date)
        out_var = resample(out_var)
        out_var = window(out_var, "latitude"; left = -60, right = 60)
        out_var = replace(out_var, NaN => NaNStatistics.nanmean(out_var.data))
        out_var = ClimaAnalysis.Var.convert_units(out_var, "kg m-2"; conversion_function = x -> x/1000)
        set_short_name!(out_var, short_name)
        out_var
    end
    return (; net_rad, sw_cre, lw_cre, rlut, rsut, rsutcs, rlutcs, rsdt, pr, ts, lwp, iwp)
end

# The ERA5 pressure range is not as large as the ClimaAtmos default pressure levels,
# so we need to limit outputvars to the ERA5 pressure range (100.0 - 100_000.0 Pa) 
limit_pressure_dim_to_era5_range(v) = limit_pressure_dim(v, 100.0, 100_000.0)

function limit_pressure_dim(output_var, min_pressure, max_pressure)
    @assert has_pressure(output_var)
    pressure_dims = output_var.dims[pressure_name(output_var)]
    valid_pressure_levels = filter(pressure_dims) do pressure
        min_pressure <= pressure <= max_pressure
    end
    lowest_valid_level = minimum(valid_pressure_levels)
    highest_valid_level = maximum(valid_pressure_levels)
    return window(output_var, "pfull"; left = lowest_valid_level, right = highest_valid_level)
end

get_monthly_averages(simdir, var_name) = get(simdir; short_name = var_name, reduction = "average", period = "1M")

"""
    seasonally_aligned_yearly_average(var, yr)

Return the yearly average of `var` for the specified year `yr`, using a seasonal 
alignment from December of the previous year to November of the given year.
"""
function seasonally_aligned_yearly_average(var, yr)
    year_window = window(var, "time"; left = DateTime(yr - 1, 12, 1), right = DateTime(yr, 11, 1))
    return average_time(year_window)
end

"""
    get_seasonal_covariance(output_var; model_error_scale = nothing, regularization = nothing)

Computes the diagonal covariance matrix of seasonal averages of `output_var`.

# Arguments
- `output_var`: Climate variable data (OutputVar or similar)
- `model_error_scale`: Optional scaling factor for model error, applied as a fraction of the mean
- `regularization`: Optional regularization term added to variance values
"""
function get_seasonal_covariance(output_var; model_error_scale = nothing, regularization = nothing)
    seasonal_averages = average_season_across_time(output_var)
    variance_per_season = map(split_by_season(seasonal_averages)) do season
        variance = flatten(variance_time(season)).data
        if !isnothing(model_error_scale)
            variance .+= (model_error_scale .* flatten(average_time(season)).data) .^ 2
        end
        return variance
    end
    diag_cov = vcat(variance_per_season...)
    !isnothing(regularization) && (diag_cov .+= regularization)
    return Diagonal(diag_cov)
end

"""
    year_of_seasonal_averages(output_var, yr)

Compute seasonal averages for a specific year `yr` from the `output_var`.
"""
function year_of_seasonal_averages(output_var, yr)
    seasonal_averages = average_season_across_time(output_var)
    season_and_years = map(ClimaAnalysis.Utils.find_season_and_year, dates(seasonal_averages))
    indices = findall(s -> s[2] == yr, season_and_years)
    isempty(indices) && error("No data found in $(long_name(output_var)) for the given year: $yr")
    min_idx, max_idx = extrema(indices)
    left = dates(seasonal_averages)[min_idx]
    right = dates(seasonal_averages)[max_idx]
    return window(seasonal_averages, "time"; left, right)
end

"""
    make_single_year_of_seasonal_observations(output_var, yr)

Create an `EKP.Observation` for a specific year `yr` from the `output_var`.
"""
function make_single_year_of_seasonal_observations(output_var, yr)
    seasonal_averages = year_of_seasonal_averages(output_var, yr)
    # Split into four OutputVars to get the same format as the covariance matrix
    obs_vec = flatten_seasonal_averages(seasonal_averages)
    obs_cov = get_seasonal_covariance(output_var; model_error_scale = 0.05, regularization = 1e-3)
    name = get(output_var.attributes, "CF_name", get(output_var.attributes, "long_name", ""))
    return EKP.Observation(obs_vec, obs_cov, "$(yr)_$name")
end

flatten_seasonal_averages(seasonal_averages) = vcat(map(x -> flatten(x).data, split_by_season(seasonal_averages))...)
# TODO: Add exception for net rad
"""
    create_observation_vector(nt, yrs = 17)

Given a NamedTuple, produce a vector of `EKP.Observation`s, where each observation
consists of seasonally averaged fields, with the exception of globally averaged yearly radiative balance
"""
function create_observation_series(nt; short_names = keys(nt), model_error_scale = 0.05, regularization = 1e-6, year_range = 2002:2018, batch_size = 1)
    nt = (;filter(p -> p.first in short_names, pairs(nt))...)
    vars = average_season_across_time.(values(nt))
    sample_dates = [
        (DateTime(i, 12, 1), DateTime(i+1, 9, 1)) for i in year_range
    ]
    cov = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(sample_dates; model_error_scale, regularization)
    obs_vec = map(sample_dates) do (start_date, end_date)
        obs = ClimaCalibrate.ObservationRecipe.observation(cov, vars, start_date, end_date)
    end
    observation_series = ClimaCalibrate.observation_series_from_samples(obs_vec,batch_size)
    return observation_series
end
