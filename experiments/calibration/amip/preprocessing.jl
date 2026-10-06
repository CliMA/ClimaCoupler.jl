# A collection of preprocessing functions
# Ideally, these should be reused for both the observational and simulation data

import ClimaAnalysis
import ClimaCoupler
import Statistics

"""
    select_pressure_levels(var, pressure_levels::Union{Vector, AbstractFloat})

Select the `pressure_levels` from `var` if pressure is a dimension of `var`.
"""
function select_pressure_levels(var, pressure_levels::Union{Vector, AbstractFloat})
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
    first_lat_idx = findfirst(lat -> lat >= lat_left, lats)
    last_lat_idx = findlast(lat -> lat <= lat_right, lats)
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
    zonal_average(var)

Average `var` over longitude, ignoring `NaN`s. A variable with no longitude dimension is
returned unchanged.

Native lon/lat points are strongly correlated, so treating them as independent
observations over-informs the inverse and collapses the ensemble. Apply this identically
to the observations (`generate_observations.jl`) and to the simulation
(`observation_map.jl`): a mismatch does not error, it silently compares one field against
a different one.
"""
function zonal_average(var)
    ClimaAnalysis.has_longitude(var) ||
        error("Variable $(ClimaAnalysis.short_name(var)) has no longitude")
    @info "Zonal (longitude) averaging $(ClimaAnalysis.short_name(var))"
    return ClimaAnalysis.average_lon(var; ignore_nan = true)
end

"""
    coverage_mask_path(output_dir)

Path of the saved observational coverage masks. Written by `generate_observations.jl` and
read by `observation_map.jl`.
"""
coverage_mask_path(output_dir) = joinpath(output_dir, "coverage_masks.jld2")

"""
    coverage_mask(var, date_ranges)

Where `var` has observational data at every date in `date_ranges`, as a
ClimaAnalysis longitude-latitude mask. Call the result on an `OutputVar` to set
its uncovered points to `NaN`. Returns `nothing` for a product that covers every
point, so a global product carries no mask at all.

Use this whenever an observation does not cover every point, and apply it to the
observation and to the simulation alike before any spatial reduction. A
reduction such as `zonal_average` ignores `NaN`, so an observation reduces over
its covered points while an unmasked simulation reduces over all of them, and
the two results are different quantities even though they have the same shape.
Applying the mask to the observation as well gives every date one fixed
coverage, so the interannual spread the covariance is estimated from is climate
rather than year-to-year coverage wobble.

This function was written to support MAC `lwp`. It is an ocean-only retrieval,
missing over land at about half of the grid points, so its zonal mean is an
ocean average while an unmasked simulation averages ocean and land together.
Modelled liquid water path over land is much lower than over ocean, so the
difference between them partly measures each band's land fraction rather than
the model's cloud.

Restricted to `date_ranges` on purpose. A record spanning decades, unioned over
every slice, drops any point that is missing in any single month.

Longitude and latitude only. A masked product with another dimension, such as a
cloud fraction on levels, errors here rather than guessing how to collapse it.
"""
function coverage_mask(var, date_ranges)
    tname = ClimaAnalysis.time_name(var)
    keep = findall(
        date -> any(range -> first(range) <= date <= last(range), date_ranges),
        ClimaAnalysis.dates(var),
    )
    isempty(keep) && error(
        "None of the date ranges $(date_ranges) are in $(ClimaAnalysis.short_name(var)); check them against the observational data.",
    )

    # One coverage for every date in use: missing at any of them is missing at all of them.
    in_use = ClimaAnalysis.select(var; by = ClimaAnalysis.Index(), time = keep)
    ClimaAnalysis.propagate_nans!(in_use)
    uncovered = isnan.(Array(selectdim(in_use.data, in_use.dim2index[tname], 1)))
    any(uncovered) || return nothing

    dims = empty(var.dims)
    dim_attributes = empty(var.dim_attributes)
    for name in filter(!=(tname), collect(keys(var.dims)))
        dims[name] = var.dims[name]
        haskey(var.dim_attributes, name) &&
            (dim_attributes[name] = var.dim_attributes[name])
    end
    @info "Coverage mask for $(ClimaAnalysis.short_name(var))" missing_fraction =
        round(count(uncovered) / length(uncovered); digits = 3)
    return ClimaAnalysis.generate_lonlat_mask(
        ClimaAnalysis.OutputVar(
            deepcopy(var.attributes),
            dims,
            dim_attributes,
            Float64.(.!uncovered),
        ),
        NaN,
        1,
    )
end

"""
    get_lonlat_regridder(config_file)

Create a regridder for `OutputVar`s for regridding to the simulation grid.
"""
function get_lonlat_regridder(config_file)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    if !isnothing(get(config_dict, "netcdf_interpolation_num_points", nothing))
        (nlon, nlat, _) = tuple(config_dict["netcdf_interpolation_num_points"]...)
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
    compute_mean_and_stddev(normalization_stas, var::ClimaAnalysis.OutputVar)

Generate normalization statistics by computing a single mean and standard
deviation for `var`, over the finite data only.
"""
function compute_mean_and_stddev(var::ClimaAnalysis.OutputVar)
    finite_data = filter(isfinite, var.data)
    isempty(finite_data) && error("$(ClimaAnalysis.short_name(var)) has no finite data")
    mean_of_var = Statistics.mean(finite_data)
    std_of_var = Statistics.std(finite_data)
    std_of_var ≈ 0.0 && error("Standard deviation is zero; check your data")
    return (mean_of_var, std_of_var)
end

"""
    compute_normalization!(normalization_stats::Dict, var)

Update `normalization_stats` with a pair of (short_name, pressure_level) to
(mean, std).

For variables without pressure levels, `pressure_level` is set to `nothing`.

Normalization statistics are computed for each variable and pressure level
combination.
"""
function compute_normalization!(normalization_stats::Dict, var)
    if ClimaAnalysis.has_pressure(var)
        for pressure_level in ClimaAnalysis.pressures(var)
            var_view_of_pressure_level = ClimaAnalysis.view_select(
                var,
                by = ClimaAnalysis.MatchValue(),
                pressure_level = pressure_level,
            )
            var_mean, var_stddev = compute_mean_and_stddev(var_view_of_pressure_level)
            normalization_stats[(
                ClimaAnalysis.short_name(var_view_of_pressure_level),
                pressure_level,
            )] = (var_mean, var_stddev)
        end
    else
        var_mean, var_stddev = compute_mean_and_stddev(var)
        normalization_stats[(ClimaAnalysis.short_name(var), nothing)] =
            (var_mean, var_stddev)
    end
    return nothing
end

"""
    apply_normalization!(normalization_stats, var::ClimaAnalysis.OutputVar)

Apply normalization using the statistics saved in `normalization_stats`.
"""
function apply_normalization!(normalization_stats, var::ClimaAnalysis.OutputVar)
    if ClimaAnalysis.has_pressure(var)
        for pressure_level in ClimaAnalysis.pressures(var)
            var_view_of_pressure_level = ClimaAnalysis.view_select(
                var,
                by = ClimaAnalysis.MatchValue(),
                pressure_level = pressure_level,
            )

            (ClimaAnalysis.short_name(var), pressure_level) in keys(normalization_stats) ||
                continue

            mean_var, std_var =
                normalization_stats[(ClimaAnalysis.short_name(var), pressure_level)]
            var_view_of_pressure_level.data .-= mean_var
            var_view_of_pressure_level.data ./= std_var
        end
    else
        (ClimaAnalysis.short_name(var), nothing) in keys(normalization_stats) || return
        mean_var, std_var = normalization_stats[(ClimaAnalysis.short_name(var), nothing)]
        var.data .-= mean_var
        var.data ./= std_var
    end
    return nothing
end
