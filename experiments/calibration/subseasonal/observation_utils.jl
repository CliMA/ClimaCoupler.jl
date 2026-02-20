import ClimaCoupler
using Statistics
import Dates
using ClimaAnalysis

# =============================================================================
# Shared variable units (used by both subseasonal and subseasonal_weekly)
# =============================================================================

var_units = Dict(
    "pr" => "kg m^-2 s^-1",
    "mslp" => "Pa",
    "tas" => "K",
    "tas - ta" => "K",
    "ts" => "K",
    "sp" => "Pa",
    "hfls" => "W m^-2",
    "hfss" => "W m^-2",
    "rsus" => "W m^-2",
    "rlus" => "W m^-2",
    "rsut" => "W m^-2",
    "rlut" => "W m^-2",
)

# =============================================================================
# Shared utility functions (used by both subseasonal and subseasonal_weekly)
# =============================================================================

"""
    remove_global_mean(var)

Subtract the latitude-weighted global mean from the given `OutputVar` `var`.
"""
function remove_global_mean(var)
    mean_var = ClimaAnalysis.average_lonlat(var; weighted = true)
    mean_data = mean_var.data
    return ClimaAnalysis.replace(val -> val - mean_data[1], var)
end

"""
    get_var(short_name, simdir)

Get an `OutputVar` from `simdir` for the given `short_name`.
Handles composite variables like "tas - ta" (surface-minus-boundary-layer temperature).
"""
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

Preprocess simulation output `var` before flattening for the G ensemble matrix.

Shifts the time axis by the largest period in the date range, sets units, and
windows to the date range. Calls `largest_period` (defined per-pipeline) to
determine the shift.
"""
function preprocess_var(var, sample_date_range)
    period = largest_period(sample_date_range)
    var = ClimaAnalysis.Var._shift_by(var, date -> date - period)
    var = set_units(var, var_units[short_name(var)])
    var = window(var, "time"; left = sample_date_range[1], right = sample_date_range[2])
    return var
end
