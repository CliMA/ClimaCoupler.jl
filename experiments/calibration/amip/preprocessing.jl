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
    select_altitude_levels(var, altitudes)

Select `altitudes` (meters) from `var` if it has an altitude dimension, using
nearest-level selection, then pin the altitude coordinate to the nominal
`altitudes`.

Cloud fraction `cl` is compared on altitude levels, but the observation
(CALIPSO, 240 m grid) and the simulation (model z-levels) are on different
vertical grids, so an exact `MatchValue` selection (as used for pressure levels)
would not line up. We instead take the nearest level on each grid and relabel the
coordinate to the shared nominal targets so obs and sim align positionally in the
flattened observation vector. The nearest levels differ by at most ~120 m between
the two grids, which is negligible relative to cloud vertical structure.

This is a no-op for variables without an altitude dimension (e.g. `lwp`).
"""
function select_altitude_levels(var, altitudes)
    ClimaAnalysis.has_altitude(var) || return var
    @info "Selecting altitude levels: $(altitudes) for $(ClimaAnalysis.short_name(var))"
    var = ClimaAnalysis.select(
        var;
        by = ClimaAnalysis.NearestValue(),
        altitude = altitudes,
    )
    # Relabel the altitude coordinate to the nominal targets so obs and sim share
    # an identical vertical axis.
    zname = ClimaAnalysis.altitude_name(var)
    new_dims = copy(var.dims)
    new_dims[zname] = collect(float.(altitudes))
    return ClimaAnalysis.remake(var; dims = new_dims)
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

Reduce `var` to a NaN-aware zonal (longitude) mean, collapsing the longitude
dimension so each latitude (and vertical level) contributes a single constraint.

Rationale: a global field's grid points are highly spatially correlated, but the
`SVDplusDCovariance` only captures a handful of interannual modes
(rank ≤ n_covariance_dates − 1) and treats every remaining grid point as an
independent, tightly-constrained observation. With O(10⁴) points that massively
over-informs a 3-parameter inverse — the effective information is ~10⁷ and the
`TransformUnscented` ensemble collapses to a point after a single update (spread
→ 1e-13), after which no further learning happens. Averaging zonally cuts the
constraint count ~2 orders of magnitude down to the true large-scale degrees of
freedom, and simultaneously averages out per-gridpoint weather noise (better
signal-to-noise). Must be applied identically to obs and sim so the flattened
observation vectors stay positionally aligned.

No-op for variables without a longitude dimension.
"""
function zonal_average(var)
    ClimaAnalysis.has_longitude(var) || return var
    @info "Zonal (longitude) averaging $(ClimaAnalysis.short_name(var))"
    return ClimaAnalysis.average_lon(var; ignore_nan = true)
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
    harmonize_nan_mask_over_dates!(var, date_ranges)

Force `var` to share a single spatial NaN pattern across all dates in
`date_ranges`: any location that is NaN at ANY of those dates is set to NaN at
every time slice.

This is required by `SVDplusDCovariance`, whose per-sample flattening drops NaNs
(`ClimaAnalysis.flatten(...; ignore_nan = true)`) and errors when the samples
differ in length. Satellite `lwp` (MAC) has coverage that varies year to year, so
without this each covariance year would drop a different number of points.
Variables with a time-invariant NaN pattern (e.g. ERA5 below-ground points) are
unaffected. The cost is that the calibration sample uses only locations with
complete coverage across all `date_ranges`.
"""
function harmonize_nan_mask_over_dates!(var, date_ranges)
    tname = ClimaAnalysis.time_name(var)
    isnothing(tname) && return var
    tidx = var.dim2index[tname]
    all_dates = ClimaAnalysis.dates(var)

    # Time indices covered by the requested date ranges (inclusive), mirroring
    # how ObservationRecipe windows each sample.
    sel = Int[]
    for (s, e) in date_ranges
        idxs = findall(d -> s <= d <= e, all_dates)
        isempty(idxs) && error(
            "Date range ($s, $e) not found in $(ClimaAnalysis.short_name(var)); " *
            "check COVARIANCE_DATE_RANGES against the observational data.",
        )
        append!(sel, idxs)
    end

    # Union of NaN locations across the selected time slices (spatial dims only).
    nanmask = falses(size(selectdim(var.data, tidx, first(sel)))...)
    for i in sel
        nanmask .|= isnan.(selectdim(var.data, tidx, i))
    end

    # Apply that union mask to every time slice.
    for i in axes(var.data, tidx)
        selectdim(var.data, tidx, i)[nanmask] .= NaN
    end
    return var
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
deviation for `var`.
"""
function compute_mean_and_stddev(var::ClimaAnalysis.OutputVar)
    # NaN-aware: variables such as `lwp` (satellite retrieval) have missing
    # points. `Statistics.mean`/`std` propagate NaN, which would make the
    # normalization constants NaN and turn the entire normalized field into NaN.
    finite_data = filter(isfinite, var.data)
    isempty(finite_data) && error("No finite data to normalize; check your data")
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
