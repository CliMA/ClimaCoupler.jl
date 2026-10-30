# A collection of preprocessing functions
# Ideally, these should be reused for both the observational and simulation data

import ClimaAnalysis
import ClimaCoupler
import Dates
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
    global_mean(var)

The area-weighted global mean of `var`, ignoring `NaN`s, as a function of time.
"""
global_mean(var) = ClimaAnalysis.average_lonlat(var; ignore_nan = true, weighted = true)

"""
    seasonal_average(var)

Average `var` over each meteorological season (DJF, MAM, JJA, SON) of each year. A season is
dated by its first month, so a calibration targets SON 2010 with the sample date range
`(DateTime(2010, 9, 1), DateTime(2010, 9, 1))` and `extend = Dates.Month(3)`. The
simulation's spinup months form their own season, which no sample date selects.
"""
function seasonal_average(var)
    @info "Seasonal averaging $(ClimaAnalysis.short_name(var))"
    return ClimaAnalysis.average_season_across_time(var)
end

"""
    observable(config, name)

The entry of `config.observables` for the loss field `name` (see `CalibrateConfig`).
"""
function observable(config, name)
    i = findfirst(o -> o.name == name, config.observables)
    isnothing(i) && error("$name is not one of the config's observables")
    return config.observables[i]
end

"""
    loss_values(vars, coverage_masks, config)

The values the loss compares, computed identically for the observations
(`generate_observations.jl`) and the model output (`observation_map.jl`). Change the loss here
and both sides change. For each field:

1. Restrict it to the points its observation covers (`coverage_masks`, see `coverage_mask`).
   A satellite retrieval is not global (MAC `lwp` is ocean only) and the averages ignore
   `NaN`, so without this the model would also average points the observation never saw.
   For `lwp` that changed the area-weighted bias by a factor of four.
2. Apply the observable's `reduce` function.
"""
function loss_values(vars, coverage_masks, config)
    return map(vars) do var
        name = ClimaAnalysis.short_name(var)
        mask = get(coverage_masks, name, nothing)
        var = isnothing(mask) ? var : mask(var)
        return observable(config, name).reduce(var)
    end
end

"""
    observation_vector_path(output_dir)

Path of the observation vector, kept with the run it belongs to. Written by
`generate_observations.jl` and read by `run_calibration.jl`.
"""
observation_vector_path(output_dir) = joinpath(output_dir, "observation_vec.jld2")

"""
    check_extend(vars, config)

Error unless members run long enough to produce one full value of every reduced field:
`extend` must cover the spacing of the reduced observations `vars`, e.g. `Month(3)` for
seasonal means. A shorter run would silently average fewer months than the observation.
"""
function check_extend(vars, config)
    start = first(first(config.sample_date_ranges))
    run_after_target = (start + config.extend) - start
    for var in vars
        dates = sort(unique(ClimaAnalysis.dates(var)))
        length(dates) < 2 && continue
        spacing = minimum(diff(dates))
        run_after_target >= spacing || error(
            "$(ClimaAnalysis.short_name(var)) has one value every \
            $(Dates.value(spacing) ÷ 86_400_000) days after its reduction, but members run \
            only $(config.extend) past each target date. Set extend to cover it, for example \
            Month(3) for seasonal means.",
        )
    end
    return nothing
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

Pass every month that goes into a reduced value, from each target date to the end of its
`extend`, not only the target dates. For SON, a point the observation misses in November
but not in September would otherwise stay unmasked. The averages ignore `NaN`, so the
observation's seasonal mean there would use September and October only while the model's
uses all three months, and a zonal mean of the season would differ from a seasonal mean of
the zonal means. With one coverage for every month, both sides average the same points in
every month, in either order.

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
