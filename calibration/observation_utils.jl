import ClimaCoupler
using Statistics
import Dates

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/setup_run.jl"))
ext = Base.get_extension(ClimaCalibrate, :ClimaAnalysisExt)

"""
    resampled_lonlat(config_file)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.
"""
function resampled_lonlat(config_file)
    cs = CoupledSimulation(config_file)
    center_space = cs.model_sims.atmos_sim.domain.center_space
    (lon_nlevels, lat_nlevels, z_nlevels) =
        ClimaDiagnostics.Writers.default_num_points(center_space)
    longitudes = range(-180, 180, lon_nlevels)
    latitudes = range(-90, 90, lat_nlevels)
    stretch = center_space.grid.vertical_grid.topology.mesh.stretch
    # TODO: Account for stretch for 3D variables and interpolate to pressure?
    z_levels = range(dz_bottom, Spaces.z_max(center_space), z_nlevels)
    return var -> resampled_to(var; lon = longitudes, lat = latitudes)
end


# TODO: Don't need the functions below since we are assuming that the
# diagnostics are producing weekly means, but we may need them in the future

"""
    compute_weekly_mean_from_daily_mean(var::OutputVar, reference_date)

Compute weekly mean from daily mean in `var`, where the weekly mean start from
`reference_date`.

If there are `NaN`s in the data, then the mean is `NaN`.
"""
function compute_weekly_mean_from_daily_mean(var::OutputVar, reference_date)
    var = ClimaAnalysis.window(var, "time", left = reference_date)
    var_dates = ClimaAnalysis.dates(var)
    reference_date in var_dates || error("$reference_date is not in $var_dates")

    # TODO: Figure out group_by

    ext.group_and_reduce_by(var, "time", group_by, reduce_by)
    return nothing
end

function compute_weekly_sum_from_daily_mean(var::OutputVar, reference_date)

    return var
end

# TODO: This function should be written in ClimaAnalysis
function shift_to_previous_day(var)
    var_dates = ClimaAnalysis.dates(var)
    var_dates .-= Dates.Day(1)
    start_date = Dates.DateTime(var.attributes["start_date"])
    shifted_times = var_dates .- start_date
    shifted_seconds = date_to_time.(start_date, shifted_times)
end
