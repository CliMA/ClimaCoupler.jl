import ClimaCoupler
using Statistics
import Dates
using ClimaAnalysis

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/setup_run.jl"))
ext = Base.get_extension(ClimaCalibrate, :ClimaAnalysisExt)


var_units = Dict(
    "pr" => "kg m^-2 s^-1",
    "mslp" => "Pa",
    "tas" => "K",
    "tas - ta" => "K",
    "hfls" => "W m^-2",
    "hfss" => "W m^-2",
    "rsus" => "W m^-2",
    "rlus" => "W m^-2",
)

function remove_global_mean(var)
    mean_var = ClimaAnalysis.Var.average_lonlat(var; weighted = true)
    mean_data = mean_var.data
    return ClimaAnalysis.replace(val -> val - mean_data[1], var)
end
