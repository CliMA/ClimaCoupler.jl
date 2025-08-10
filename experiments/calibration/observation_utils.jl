import ClimaCoupler
using Statistics
import Dates
using ClimaAnalysis

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/setup_run.jl"))
ext = Base.get_extension(ClimaCalibrate, :ClimaAnalysisExt)

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

        function group_by(dim)
        grouped_dates =
            ClimaAnalysis.Utils.time_to_date.(Dates.DateTime(2008), dim) |>
            ClimaAnalysis.Utils.split_by_season_across_time
        grouped_times = [
            ClimaAnalysis.Utils.date_to_time.(Dates.DateTime(2008), x) for
            x in grouped_dates
        ]
        return grouped_times
    end

    ext.group_and_reduce_by(var, "time", group_by, reduce_by)
    return nothing
end