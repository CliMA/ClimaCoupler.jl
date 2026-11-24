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

The area-weighted global mean of `var`, ignoring `NaN`s. For an `OutputVar` on time,
latitude, and longitude, this is a function of time.
"""
global_mean(var) = ClimaAnalysis.average_lonlat(var; ignore_nan = true, weighted = true)

"""
    remove_global_mean(var)

`var` minus its area-weighted global mean at each time, ignoring `NaN`s: the spatial pattern
of the field without its global offset. Compose it after the time reduction, for example
`remove_global_mean ∘ seasonal_average` for each season's map minus that season's global mean.
"""
function remove_global_mean(var)
    mean = global_mean(var)
    data = copy(var.data)
    position = var.dim2index[ClimaAnalysis.time_name(var)]
    for (t, slice) in enumerate(eachslice(data; dims = position))
        slice .-= mean.data[t]
    end
    return ClimaAnalysis.OutputVar(var.attributes, var.dims, var.dim_attributes, data)
end

"""
    seasonal_average(var)

Average `var` over each meteorological season (DJF, MAM, JJA, SON) of each year. A season is
dated by its first month, and `snap_date_ranges` shrinks each `(start, stop)` target range
to the seasonal dates inside it, so `(DateTime(2009, 12, 1), DateTime(2010, 11, 1))` with
`extend = Dates.Month(1)` selects all four seasons of that December-to-November year, and
`(DateTime(2010, 9, 1), DateTime(2010, 11, 1))` with `extend = Dates.Month(1)` (or
`(DateTime(2010, 9, 1), DateTime(2010, 9, 1))` with `extend = Dates.Month(3)`) selects SON
2010. The simulation's spinup months form their own season, which no sample date selects.
"""
function seasonal_average(var)
    @info "Seasonal averaging $(ClimaAnalysis.short_name(var))"
    return ClimaAnalysis.average_season_across_time(var)
end

function complete_season_years(seasons)
    dates = ClimaAnalysis.dates(seasons)
    year_of(date) = Dates.year(date + Dates.Month(1))
    years = [findall(d -> year_of(d) == y, dates) for y in unique(year_of.(dates))]
    kept = filter(indices -> length(indices) == 4, years)
    isempty(kept) &&
        error("$(ClimaAnalysis.short_name(seasons)) has no year with all four seasons")
    return kept
end

"""
    seasonal_anomaly(var)

Each season's mean minus its year's mean, at every point: the seasonal cycle without the
annual-mean state, so a loss on it ignores the time-mean bias. A year runs from December to
November, so its four seasons are DJF, MAM, JJA, and SON, and the year's mean is the mean of
the four (`annual_average`). Years without all four seasons are dropped. The four anomalies of
a year sum to zero.
"""
function seasonal_anomaly(var)
    seasons = seasonal_average(var)
    position = seasons.dim2index[ClimaAnalysis.time_name(seasons)]
    years = complete_season_years(seasons)
    data = copy(seasons.data)
    for indices in years
        mean = sum(selectdim(seasons.data, position, k) for k in indices) ./ 4
        for k in indices
            selectdim(data, position, k) .-= mean
        end
    end
    anomalies = ClimaAnalysis.OutputVar(
        seasons.attributes,
        seasons.dims,
        seasons.dim_attributes,
        data,
    )
    return ClimaAnalysis.select(
        anomalies;
        by = ClimaAnalysis.Index(),
        time = reduce(vcat, years),
    )
end

"""
    annual_average(var)

The mean of each December-to-November year's four seasonal means (DJF, MAM, JJA, SON), dated
by its December, at every point. Years without all four seasons are dropped.

With `annual_global_mean`, `annual_pattern`, and `seasonal_anomaly`, it splits a year of
seasonal maps into three parts that add up to the maps and do not overlap, so a loss on all
three cannot trade one against another.
"""
function annual_average(var)
    seasons = seasonal_average(var)
    position = seasons.dim2index[ClimaAnalysis.time_name(seasons)]
    years = complete_season_years(seasons)
    # The first season of each year is its DJF, so the year keeps that slice's date.
    # Slices `k:k` keep the time dimension, so `cat` stacks the years along it.
    means = [
        sum(selectdim(seasons.data, position, k:k) for k in indices) ./ 4 for
        indices in years
    ]
    kept = ClimaAnalysis.select(seasons; by = ClimaAnalysis.Index(), time = first.(years))
    return ClimaAnalysis.OutputVar(
        kept.attributes,
        kept.dims,
        kept.dim_attributes,
        cat(means...; dims = position),
    )
end

"""
    annual_global_mean(var)

The area-weighted global annual mean of `var` for each December-to-November year.
"""
annual_global_mean(var) = global_mean(annual_average(var))

"""
    annual_pattern(var)

Each December-to-November year's annual map minus its area-weighted global mean.
"""
annual_pattern(var) = remove_global_mean(annual_average(var))

"""
    mean_pattern_seasonal_observables(
        name;
        field = name,
        global_mean_sd,
        annual_pattern_fraction,
        seasonal_cycle_fraction,
        effective_entries = 29,
    )

The three non-overlapping observables of `name` over a December-to-November year: its
global annual mean (`annual_global_mean`), its annual spatial pattern (`annual_pattern`),
and its four seasonal anomalies (`seasonal_anomaly`). `global_mean_sd` is the per-entry
error scale of the global mean in the field's physical units, and
`annual_pattern_fraction` and `seasonal_cycle_fraction` are the model error of the two maps
as fractions of each reduced map's time mean across the covariance dates (with the squared
mean floored at its spatial median).

`effective_entries` is one count for each of the three parts (default `29`), a 3-tuple of
counts `(n_global, n_annual, n_seasonal)`, or `nothing`. Scaling each part's covariance by
`entries / effective_entries` first divides by `entries` (normalizing the 1-entry global
mean, the 1-map annual pattern, and the 4-map seasonal cycle to their mean squared whitened
residual per entry) and then weights each part by `effective_entries` in both the loss and
`EffectiveEntriesController`. Giving all three parts `29` (from the ~175 effective spatial
entries of a split-half calibration test across the six maps of a 3-field run, and because
the four seasonal anomalies at each point sum to zero and are dominated by the annual
harmonic) gives the global mean, annual pattern, and seasonal cycle equal 1/3 block weight
without rescaling `global_mean_sd`.
"""
function mean_pattern_seasonal_observables(
    name;
    field = name,
    global_mean_sd,
    annual_pattern_fraction,
    seasonal_cycle_fraction,
    effective_entries = 29,
)
    (n_global, n_annual, n_seasonal) =
        effective_entries isa NTuple{3, Any} ? effective_entries :
        (effective_entries, effective_entries, effective_entries)
    return [
        ClimaCoupler.CalibrationTools.Observable(;
            name = "$(name)_global_mean",
            field,
            reduce = annual_global_mean,
            noise_sd = global_mean_sd,
            effective_entries = n_global,
        ),
        ClimaCoupler.CalibrationTools.Observable(;
            name = "$(name)_annual_pattern",
            field,
            reduce = annual_pattern,
            model_error_fraction = annual_pattern_fraction,
            effective_entries = n_annual,
        ),
        ClimaCoupler.CalibrationTools.Observable(;
            name = "$(name)_seasonal_cycle",
            field,
            reduce = seasonal_anomaly,
            model_error_fraction = seasonal_cycle_fraction,
            effective_entries = n_seasonal,
        ),
    ]
end

"""
    renamed(var, name)

`var` with its short name set to `name`.
"""
function renamed(var, name)
    attributes = copy(var.attributes)
    attributes["short_name"] = name
    return ClimaAnalysis.OutputVar(attributes, var.dims, var.dim_attributes, var.data)
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
    mask_and_reduce(vars, coverage_masks, config; names = nothing)

Mask each field to its observation's coverage and apply each observable's `reduce`. The
observations (`generate_observations.jl`) and the model output (`observation_map.jl`) both
go through this function, so the loss compares the same quantities on both sides. `vars` are
fields, and the result has one value per matching `(observable, var)` pair in the order of
`config.observables` (or `names`):

1. Restrict the field to the points its observation covers (`coverage_masks`, see
   `coverage_mask`). A satellite retrieval is not always global (MAC `lwp` is ocean only)
   and the averages ignore `NaN`, so without this the model would also average points the
   observation never saw.
2. Apply the observable's `reduce` function, and give the result the observable's name.
"""
function mask_and_reduce(vars, coverage_masks, config; names = nothing)
    selected =
        isnothing(names) ? config.observables :
        [observable(config, name) for name in names]
    values = ClimaAnalysis.OutputVar[]
    for o in selected
        matched = false
        for var in vars
            field = ClimaAnalysis.short_name(var)
            field == o.field || continue
            matched = true
            mask = get(coverage_masks, field, nothing)
            reduced = o.reduce(isnothing(mask) ? var : mask(var))
            push!(values, field == o.name ? reduced : renamed(reduced, o.name))
        end
        matched || error("No variable for field $(o.field) (observable $(o.name))")
    end
    return values
end

"""
    observation_vector_path(output_dir)

Path of the observation vector, kept with the run it belongs to. Written by
`generate_observations.jl` and read by `run_calibration.jl`.
"""
observation_vector_path(output_dir) = joinpath(output_dir, "observation_vec.jld2")

"""
    snap_date_ranges(date_ranges, var)

Each `(start, stop)` range in `date_ranges` shrunk to the first and last date `var` has
inside `[start, stop]`. `SampleBuilder.build_samples_by_times` needs both endpoints of a
range to match exact timestamps in `var`, and a time reduction leaves fewer timestamps than
the monthly input window spans: `annual_average` over December 2009 to November 2010 is
dated by its first slice `2009-12-01`, and `seasonal_anomaly` has four seasonal timestamps
`2009-12-01` through `2010-09-01`.

Errors if any range has no timestamps in `var` or if different ranges in `date_ranges`
contain different numbers of timestamps (for example when a range is only partially covered
by the observational record).
"""
function snap_date_ranges(date_ranges, var)
    dates = sort(unique(ClimaAnalysis.dates(var)))
    counts = Int[]
    snapped = map(date_ranges) do (start, stop)
        inside = filter(date -> start <= date <= stop, dates)
        isempty(inside) &&
            error("$(ClimaAnalysis.short_name(var)) has no values between $start and $stop")
        push!(counts, length(inside))
        (first(inside), last(inside))
    end
    allequal(counts) || error(
        "$(ClimaAnalysis.short_name(var)) has $counts timestamps across $date_ranges, expected the same count in every range",
    )
    return snapped
end

"""
    check_extend(vars, config)

Error unless, for every target range `(start, stop)` in `config.sample_date_ranges`,
members run long enough to produce the reduced values in `vars` and `(start, stop)` selects
every reduced value that the simulated target window `[start, stop + extend)` produces. A
shorter `extend` would average fewer months than the observation, while a single-date range
`(start, start)` paired with a multi-period `extend` (such as `Month(12)` on
`seasonal_anomaly`) would simulate four seasons and silently select only the first.
"""
function check_extend(vars, config)
    for (start, stop) in unique(config.sample_date_ranges)
        window_end = stop + config.extend - Dates.Day(1)
        run_after_target = (stop + config.extend) - start
        for var in vars
            dates = sort(unique(ClimaAnalysis.dates(var)))
            if length(dates) >= 2
                spacing = minimum(diff(dates))
                run_after_target >= spacing || error(
                    "$(ClimaAnalysis.short_name(var)) has one value every \
                    $(Dates.value(spacing) ÷ 86_400_000) days after its reduction, but members run \
                    only $(config.extend) past each target date. Set extend to cover it, for example \
                    Month(3) for seasonal means.",
                )
            end
            selected = filter(d -> start <= d <= stop, dates)
            simulated = filter(d -> start <= d <= window_end, dates)
            if !isempty(simulated) && length(selected) != length(simulated)
                error(
                    "$(ClimaAnalysis.short_name(var)) has $(length(simulated)) value(s) in the " *
                    "simulated target window ($start to $(stop + config.extend)), but " *
                    "sample_date_ranges ($start, $stop) selects $(length(selected)) of them. " *
                    "Set the range end date to $(last(simulated)) and extend to one averaging period.",
                )
            end
        end
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
        "None of the date ranges $(date_ranges) are in $(ClimaAnalysis.short_name(var)), check them against the observational data.",
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
