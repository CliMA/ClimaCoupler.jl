import ClimaCoupler
using Statistics
import Dates
using ClimaAnalysis

# Shared variable units (used by both subseasonal and subseasonal_weekly)

var_units = Dict(
    # Surface temperature/pressure
    "pr" => "kg m^-2 s^-1",
    "mslp" => "Pa",
    "tas" => "K",
    "tas - ta" => "K",
    "ts" => "K",
    "sp" => "Pa",
    # Surface fluxes
    "hfls" => "W m^-2",
    "hfss" => "W m^-2",
    # Radiation (TOA)
    "rsdt" => "W m^-2",    # incoming solar
    "rsut" => "W m^-2",    # TOA outgoing SW
    "rlut" => "W m^-2",    # TOA outgoing LW
    "rsutcs" => "W m^-2",  # TOA outgoing SW clear-sky
    "rlutcs" => "W m^-2",  # TOA outgoing LW clear-sky
    # Radiation (surface)
    "rsds" => "W m^-2",    # surface downwelling SW
    "rsus" => "W m^-2",    # surface upwelling SW
    "rlds" => "W m^-2",    # surface downwelling LW
    "rlus" => "W m^-2",    # surface upwelling LW
    "rsdscs" => "W m^-2",  # surface downwelling SW clear-sky
    "rsuscs" => "W m^-2",  # surface upwelling SW clear-sky
    "rldscs" => "W m^-2",  # surface downwelling LW clear-sky
    "ta" => "K",      # Temperature at any pressure level
    "hur" => "unitless",  # Relative humidity at any pressure level
    "hus" => "unitless",  # Specific humidity at any pressure level
)


const PRESSURE_LEVEL_BASE_UNITS = Dict(
    "ta" => "K",
    "hur" => "unitless",
    "hus" => "unitless",
)

"""
    get_var_units(short_name)

Get the units for a variable, handling pressure-level variables automatically.
For "ta_850hPa", returns units for "ta", etc.
"""
function get_var_units(short_name::String)

    if haskey(var_units, short_name)
        return var_units[short_name]
    end

    if is_pressure_level_variable(short_name)
        base_name, _ = parse_pressure_level_variable(short_name)
        if haskey(PRESSURE_LEVEL_BASE_UNITS, base_name)
            return PRESSURE_LEVEL_BASE_UNITS[base_name]
        end
    end
    
    error("Unknown variable: $short_name - add it to var_units or PRESSURE_LEVEL_BASE_UNITS")
end

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
    is_pressure_level_variable(short_name)

Check if a short_name represents a pressure-level variable (e.g., "ta_850hPa").
"""
function is_pressure_level_variable(short_name::String)
    return occursin(r"_\d+hPa$", short_name)
end

"""
    parse_pressure_level_variable(short_name)

Parse a pressure-level variable name into (base_name, pressure_level).
e.g., "ta_850hPa" -> ("ta", 850.0)
"""
function parse_pressure_level_variable(short_name::String)
    m = match(r"^(.+)_(\d+)hPa$", short_name)
    if isnothing(m)
        error("Invalid pressure-level variable name: $short_name")
    end

    return (String(m.captures[1]), parse(Float64, m.captures[2]))
end

"""
    slice_pressure_level(var, pressure_hPa)

Slice an `OutputVar` at the given pressure level (in hPa), converting to Pa
if the data's pressure dimension is stored in Pa.

Uses ClimaAnalysis's built-in `pressure_name` to resolve the dimension
(handles both "pfull" and "pressure_level"), and `slice` already aliases
between these names, so we just need to pass the correct numeric value.
"""
function slice_pressure_level(var, pressure_hPa)
    pdim = ClimaAnalysis.pressure_name(var)
    dim_unit = get(get(var.dim_attributes, pdim, Dict()), "units", "")

    if dim_unit == "Pa"
        pressure_val = pressure_hPa * 100.0
    elseif dim_unit == "hPa"
        pressure_val = pressure_hPa
    else
        error("Unknown pressure units '$dim_unit' for dimension '$pdim'. " *
              "Expected 'Pa' or 'hPa'.")
    end

    return ClimaAnalysis.slice(var, pfull = pressure_val, by = ClimaAnalysis.MatchValue())
end

"""
    get_var(short_name, simdir)

Get an `OutputVar` from `simdir` for the given `short_name`.
Handles:
- Composite variables like "tas - ta" (surface-minus-boundary-layer temperature)
- Pressure-level variables like "ta_850hPa" (temperature at 850 hPa)
"""
function get_var(short_name, simdir)
    if short_name == "tas - ta"
        tas = get(simdir; short_name = "tas")
        ta = get(simdir; short_name = "ta")
        ta_900hpa = slice(ta; z = 1000)
        var = tas - ta_900hpa
    elseif is_pressure_level_variable(short_name)

        base_name, pressure_hPa = parse_pressure_level_variable(short_name)

        var = get(simdir; short_name = base_name, coord_type = "pressure")
        var = slice_pressure_level(var, pressure_hPa)

        # Handle relative humidity units
        # ClimaDiagnostics outputs hur with empty string units (""), but
        # ClimaAnalysis treats "" as missing units and expects "unitless".
        if base_name == "hur"
            units = ClimaAnalysis.units(var)
            if units == "%"
                var = ClimaAnalysis.convert_units(
                    var, "unitless";
                    conversion_function = x -> 0.01 * x,
                )
            elseif units == "" || units == "1"
                var = ClimaAnalysis.set_units(var, "unitless")
            end
        end
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
    var = set_units(var, get_var_units(short_name(var)))
    var = window(var, "time"; left = sample_date_range[1], right = sample_date_range[2])
    return var
end
