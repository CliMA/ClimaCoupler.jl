# Utilities for processing observations/OutputVars
# Used to generate observations and compute the observation map
using ClimaAnalysis, ClimaCoupler
push!(ClimaAnalysis.Var.ALTITUDE_NAMES, "height")
using OrderedCollections
import ClimaUtilities.ClimaArtifacts: @clima_artifact
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
const FT = Float32

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/data_sources.jl"))



"""
    create_observation_vector(nt, yrs = 17)

Given a NamedTuple, produce a vector of `EKP.Observation`s, where each
observation consists of seasonally averaged fields, with the exception of
globally averaged yearly radiative balance
"""
function create_observation_series(
    nt;
    short_names = keys(nt),
    model_error_scale = 0.05,
    regularization = 1e-6,
    year_range = 2006:2017,
    batch_size = 1,
)
    net_rad = nt.net_rad
    net_rad_var = [variance_time(net_rad).data...;;]
    nt = (; filter(p -> p.first in short_names, pairs(nt))...)
    vars = average_season_across_time.(values(nt))
    # TODO: Change this to monthly averages
    sample_dates = [(DateTime(i, 12, 1), DateTime(i + 1, 9, 1)) for i in year_range]
    cov = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(sample_dates; model_error_scale, regularization)
    obs_vec = map(sample_dates) do (start_date, end_date)
        @info "Creating observations for $(start_date, end_date)"
        # net_rad_year = window(net_rad, "time"; left = start_date, right = start_date + Month(11))
        # @assert length(dates(net_rad_year)) == 12
        # net_rad_obs = EKP.Observation(
        #     ClimaAnalysis.flatten(average_time(net_rad_year)).data,
        #     net_rad_var,
        #     "net_rad_$(string(start_date))",
        # )
        ClimaCalibrate.ObservationRecipe.observation(cov, vars, start_date, end_date)
        # obs = EKP.combine_observations([net_rad_obs, obs])
    end
    observation_series = ClimaCalibrate.observation_series_from_samples(obs_vec, batch_size)
    return observation_series
end

"""
    get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)

Return a NamedTuple of OutputVars containing all initial coarse AMIP
observations. Start date is set to `DateTime(2000, 3, 1)`. All OutputVars are
resampled to the model diagnostic grid.
"""
function get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)
    if !isnothing(diagnostic_var3d)
        diagnostic_var3d = limit_pressure_to_era5_range(diagnostic_var3d)
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

    # # TOA net radiative flux
    net_rad = rlut + rsut - rsdt |> average_lonlat
    set_short_name!(net_rad, "net_rad")
    # cloud radiative effect
    sw_cre = rsut - rsutcs
    set_short_name!(sw_cre, "sw_cre")
    lw_cre = rlut - rlutcs
    set_short_name!(lw_cre, "lw_cre")
    # pr = resample(rad_and_pr_obs_dict["pr"](start_date))
    # set_short_name!(pr, "pr")
    # ta = resample(era5_outputvar(joinpath(obs_dir, "era5_monthly_avg_pressure_level_t.nc")))
    # set_short_name!(ta, "ta")

    # ts = resample(
    #     OutputVar(
    #         joinpath(obs_dir, "era5_monthly_avg_ts.nc");
    #         new_start_date = start_date,
    #         shift_by = Dates.firstdayofmonth,
    #     ),
    # )
    # set_short_name!(ts, "ts")

    lwp, iwp = map(["lwp", "iwp"]) do short_name
        out_var = OutputVar(joinpath(pkgdir(ClimaCoupler), "modis_lwp_iwp.nc"), short_name, new_start_date = start_date)
        out_var = resample(out_var)
        out_var = window(out_var, "latitude"; left = -60, right = 60)
        out_var = replace(out_var, NaN => NaNStatistics.nanmean(out_var.data))
        out_var = ClimaAnalysis.Var.convert_units(out_var, "kg m-2"; conversion_function = x -> x / 1000)
        set_short_name!(out_var, short_name)
        out_var
    end

    # cl = OutputVar(
    #     joinpath(@clima_artifact("calipso_cloudsat"), "radarlidar_seasonal_2.5x2.5.nc"),
    #     "cloud_fraction_on_levels",
    #     new_start_date = start_date,
    # )
    # cl = replace(cl, missing => NaN)
    # @assert size(cl.data)[4] == 2
    # new_dims = OrderedDict(k => v for (k, v) in cl.dims if k != "doop")
    # new_dim_attributes = OrderedDict(k => v for (k, v) in cl.dim_attributes if k != "doop")
    # # TODO: use the `by` keyword argument for `slice` using `ClimaAnalysis.Index()`to slice out the doop dimension
    # cl = remake(cl; data = cl.data[:, :, :, 1, :], dims = new_dims, dim_attributes = new_dim_attributes)
    # cl = resample(cl)
    # set_short_name!(cl, "cl")

    # # TODO: Ensure latlon grid is aligned with coupler output
    # pr_bin_data = JLD2.load_object("precip_seasonal_histograms.jld2")

    # # Split bins apart
    # concatenated_bins = Dict{Int, Array{Int64, 3}}()
    # for bin in 1:5
    #     # Collect all data for this bin across all years and seasons
    #     bin_data = []

    #     for year in keys(pr_bin_data)
    #         for season in keys(pr_bin_data[year])
    #             push!(bin_data, pr_bin_data[year][season][:, :, bin])
    #         end
    #     end
    #     concatenated_bins[bin] = cat(bin_data..., dims = 3)
    # end
    # # Get dates
    # bin_dates = Date[]
    # for year in keys(pr_bin_data)
    #     for season in keys(pr_bin_data[year])
    #         ssn_to_mnth = Dict("MAM" => 3, "JJA" => 6, "SON" => 9)
    #         if season == "DJF"
    #             push!(bin_dates, Date(year, 12, 1))  # December of the year
    #         else
    #             push!(bin_dates, Date(year + 1, ssn_to_mnth[season], 1))  # Following year for other seasons
    #         end
    #     end
    # end
    # sort!(bin_dates)
    # bin_times = Float64[Dates.Second(d - Date(start_date)).value for d in bin_dates]

    # pr_bin_dims = OrderedDict("lat" => FT.(-90.0:90), "lon" => FT.(0.0:359), "time" => bin_times)

    # lon_dict = Dict("units" => "degrees_east")
    # lat_dict = Dict("units" => "degrees_north")
    # time_dict = Dict("units" => "s")
    # pr_bin_dim_attibs = OrderedDict("lon" => lon_dict, "lat" => lat_dict, "time" => time_dict)
    # str_start_date = string(start_date)
    # #! format: off
    # pr_bin_1 = resample(OutputVar(Dict("short_name" => "pr_binned_bin_1", "start_date" =>str_start_date), pr_bin_dims, pr_bin_dim_attibs, concatenated_bins[1]))
    # pr_bin_2 = resample(OutputVar(Dict("short_name" => "pr_binned_bin_2", "start_date" =>str_start_date), pr_bin_dims, pr_bin_dim_attibs, concatenated_bins[2]))
    # pr_bin_3 = resample(OutputVar(Dict("short_name" => "pr_binned_bin_3", "start_date" =>str_start_date), pr_bin_dims, pr_bin_dim_attibs, concatenated_bins[3]))
    # pr_bin_4 = resample(OutputVar(Dict("short_name" => "pr_binned_bin_4", "start_date" =>str_start_date), pr_bin_dims, pr_bin_dim_attibs, concatenated_bins[4]))
    # pr_bin_5 = resample(OutputVar(Dict("short_name" => "pr_binned_bin_5", "start_date" =>str_start_date), pr_bin_dims, pr_bin_dim_attibs, concatenated_bins[5]))
    return (; net_rad, sw_cre, lw_cre, rlut, rsut, rsutcs, rlutcs, lwp, iwp)
    #! format: on
end

# The ERA5 pressure range is not as large as the ClimaAtmos default pressure levels,
# so we need to limit outputvars to the ERA5 pressure range (100.0 - 100_000.0 Pa) 
limit_pressure_to_era5_range(v) = limit_pressure_range(v, 100.0, 100_000.0)

function limit_pressure_range(output_var, min_pressure, max_pressure)
    @assert has_pressure(output_var)
    pressure_dims = output_var.dims[pressure_name(output_var)]
    valid_pressure_levels = filter(pressure_dims) do pressure
        min_pressure <= pressure <= max_pressure
    end
    lowest_valid_level = minimum(valid_pressure_levels)
    highest_valid_level = maximum(valid_pressure_levels)
    return window(output_var, pressure_name(output_var); left = lowest_valid_level, right = highest_valid_level)
end

function limit_altitude_range(output_var, min_level = 100.0)
    @assert has_altitude(output_var)
    z_name = altitude_name(output_var)
    altitude_dims = output_var.dims[z_name]

    valid_pressure_levels = filter(a -> min_level <= a, altitude_dims)

    lowest_valid_level = minimum(valid_pressure_levels)
    return window(output_var, z_name; left = lowest_valid_level)
end

get_monthly_averages(simdir, var_name) = get(simdir; short_name = var_name, reduction = "average", period = "1M")

set_short_name!(var, short_name) = (var.attributes["short_name"] = short_name)


weekly = OutputVar("weekly/tas_sevenday_WEEKLYMEAN_2005.nc"; new_start_date = start_date)
weekly.attributes["start_date"] = string(start_date)