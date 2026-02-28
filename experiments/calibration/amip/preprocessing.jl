# A collection of preprocessing functions
# Ideally, these should be reused for both the observational and simulation data

import ClimaAnalysis
import ClimaCoupler
import Statistics

"""
    select_pressure_levels(var, pressure_levels::Vector)

Select the `pressure_levels` from `var` if pressure is a dimension of `var`.
"""
function select_pressure_levels(var, pressure_levels::Vector)
    if ClimaAnalysis.has_pressure(var)
        @info "Selecting pressure levels: $(pressure_levels) for $(ClimaAnalysis.short_name(var))"
        var = ClimaAnalysis.select(
            var;
            by = ClimaAnalysis.MatchValue(),
            pressure_level = pressure_levels,
        )
    end
    return var
end

"""
    apply_lat_window(var, lat_left, lat_right)

Apply latitude window by constraining the longitudes to be in the range
[lat_left, lat_right].
"""
function apply_lat_window(var, lat_left, lat_right)
    lats = ClimaAnalysis.latitudes(var)
    first_lat_idx = findfirst(lon -> lon >= lat_left, lats)
    last_lat_idx = findlast(lon -> lon <= lat_right, lats)
    var = ClimaAnalysis.window(
        var,
        "latitude",
        by = ClimaAnalysis.Index(),
        left = first_lat_idx,
        right = last_lat_idx,
    )
    @info "Windowing latitudes, latitudes of $(ClimaAnalysis.short_name(var)) is $(ClimaAnalysis.latitudes(var))"
    return var
end

"""
    get_lonlat_regridder(config_file)

Create a regridder for `OutputVar`s for regridding to the simulation grid.
"""
function get_lonlat_regridder(config_file)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    if !isnothing(get(config_dict, "netcdf_interpolation_num_points", nothing))
        (nlon, nlat, nlev) = tuple(config_dict["netcdf_interpolation_num_points"]...)
    else
        # Compute from h_elem (spectral element grid)
        h_elem = get(config_dict, "h_elem", 12)
        # Default formula: h_elem * 4 panels * 3 (spectral degree)
        nlon = h_elem * 4 * 3
        nlat = nlon ÷ 2
        @info "Using model grid from h_elem=$h_elem: $(nlon)×$(nlat)"
    end
    lon_vals = range(-180, 180, nlon)
    lat_vals = range(-90, 90, nlat)
    return var -> ClimaAnalysis.resampled_as(var; longitude = lon_vals, latitude = lat_vals)
end

"""
    set_unitless_units!(var)

Set the units of `var` to "unitless" if the units is the empty string.
"""
function set_unitless_units!(var)
    if ClimaAnalysis.units(var) == ""
        # TODO: In ClimaAnalysis, there should be a set_units! function
        var.attributes["units"] = "unitless"
    end
    return var
end

"""
    compute_mean_and_and_stddev(normalization_stas, var::ClimaAnalysis.OutputVar)

Generate normalization statistics by computing a single mean and standard
deviation for `var`.
"""
function compute_mean_and_and_stddev(var::ClimaAnalysis.OutputVar)
    mean_of_var = Statistics.mean(var.data)
    std_of_var = Statistics.std(var.data)
    std_of_var ≈ 0.0 && error("Standard deviation is 0.0; check your data")
    return (mean_of_var, std_of_var)
end

"""
    compute_normalization_stats!(normalization_stats::Dict, var)

Return a dictionary mapping (short_name, pressure_level) to (mean, std).

If there is no pressure level, then the `pressure_level` is nothing.
"""
function compute_normalization_stats!(normalization_stats::Dict, var)
    if ClimaAnalysis.has_pressure(var)
        for pressure_level in ClimaAnalysis.pressures(var)
            var_view_of_pressure_level = view_select(
                var,
                by = ClimaAnalysis.MatchValue(),
                pressure_level = pressure_level,
            )
            var_mean, var_stddev = compute_mean_and_and_stddev(var)
            normalization_stats[(ClimaAnalysis.short_name(var), pressure_level)] =
                (var_mean, var_stddev)
        end
    else
        var_mean, var_stddev = compute_mean_and_and_stddev(var)
        normalization_stats[(ClimaAnalysis.short_name(var), nothing)] =
            (var_mean, var_stddev)
    end
    return nothing
end

function apply_normalization_stats!(var::ClimaAnalysis.OutputVar, normalization_stats)
    if ClimaAnalysis.has_pressure(var)
        for pressure_level in ClimaAnalysis.pressures(var)
            var_view_of_pressure_level = view_select(
                var,
                by = ClimaAnalysis.MatchValue(),
                pressure_level = pressure_level,
            )
            mean_var, std_var =
                normalization_stats[(ClimaAnalysis.short_name(var), pressure_level)]
            var_view_of_pressure_level.data .-= mean_var
            var_view_of_pressure_level.data ./= std_var
        end
    else
        mean_var, std_var = normalization_stats[(ClimaAnalysis.short_name(var), nothing)]
        var.data .-= mean_var
        var.data ./= std_var
    end
    return nothing
end
