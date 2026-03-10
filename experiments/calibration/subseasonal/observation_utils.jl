import ClimaCoupler
using Statistics
import Dates
using ClimaAnalysis

include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "code_loading.jl"))

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

"""
    remove_global_mean(var)

Subtract the latitude-weighted global mean from the given `OutputVar` `var`.
"""
function remove_global_mean(var)
    mean_var = ClimaAnalysis.average_lonlat(var; weighted = true)
    mean_data = mean_var.data
    return ClimaAnalysis.replace(val -> val - mean_data[1], var)
end
