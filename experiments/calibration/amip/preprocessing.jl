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
    lat_window()

The latitude range to keep, `(lat_left, lat_right)`. A config that defines
`LAT_WINDOW` restricts the calibration to that band, for example
`(-65, -45)` for the Southern Ocean. The whole globe otherwise.

Both the observation and the simulation must use this, or their flattened
vectors stop lining up.
"""
lat_window() = isdefined(Main, :LAT_WINDOW) ? Main.LAT_WINDOW : (-90, 90)

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
    coverage_mask_path(output_dir)

Path of the saved observational coverage mask (see `coverage_mask` /
`apply_coverage_mask`).
"""
coverage_mask_path(output_dir) = joinpath(output_dir, "coverage_masks.jld2")

"""
    ordered_dim_names(var)

Dimension names of `var` in the order of its data axes.
"""
function ordered_dim_names(var)
    idx_name = sort([(idx, name) for (name, idx) in var.dim2index])
    return [name for (_, name) in idx_name]
end

"""
    canonical_dim_name(name)

Map the various spellings of the spatial dimensions onto canonical names, so an
observation's coverage mask can be matched against a simulation variable even when
the two datasets name the same axis differently (e.g. the CALIPSO product calls the
altitude axis `height` while the model diagnostics call it `z` — without this the
mask was silently skipped with a "dimensions do not match" warning).
"""
function canonical_dim_name(name)
    n = lowercase(String(name))
    n in ("z", "height", "altitude", "z_reference") && return "altitude"
    n in ("lon", "long", "longitude") && return "longitude"
    n in ("lat", "latitude") && return "latitude"
    return n
end

"""
    date_range_nan_mask(var, date_ranges)

Union of `var`'s missing (NaN) points over the time slices covered by
`date_ranges`: `true` where the value is missing at ANY of those dates.

Restricted to `date_ranges` on purpose. Unioning over *every* slice in the record
would be wrong: satellite products span decades, so a point missing in any single
month would be dropped and the mask saturates to "everything missing".
"""
function date_range_nan_mask(var, date_ranges)
    tname = ClimaAnalysis.time_name(var)
    tidx = var.dim2index[tname]
    all_dates = ClimaAnalysis.dates(var)

    sel = Int[]
    for (s, e) in date_ranges
        idxs = findall(d -> s <= d <= e, all_dates)
        isempty(idxs) && error(
            "Date range ($s, $e) not found in $(ClimaAnalysis.short_name(var)); " *
            "check the date ranges against the observational data.",
        )
        append!(sel, idxs)
    end

    nanmask = falses(size(selectdim(var.data, tidx, first(sel)))...)
    for i in sel
        nanmask .|= isnan.(selectdim(var.data, tidx, i))
    end
    return nanmask
end

"""
    coverage_mask(var, date_ranges)

Return `(dim_names, mask)` describing where `var` has NO observational data over
`date_ranges`: `mask` is a `Bool` array over `var`'s non-time dimensions
(true = missing), and `dim_names` are those dimensions in the mask's axis order.

Used to make the simulation sample exactly the same points as the observation
before taking a zonal mean — see `apply_coverage_mask`. Shares
`date_range_nan_mask` with `harmonize_nan_mask_over_dates!`, so the mask saved for
the simulation is by construction the same one applied to the observation.
"""
function coverage_mask(var, date_ranges)
    tname = ClimaAnalysis.time_name(var)
    names = ordered_dim_names(var)
    if isnothing(tname) || !(tname in names)
        return (names, isnan.(Array(var.data)))
    end
    return (filter(!=(tname), names), date_range_nan_mask(var, date_ranges))
end

"""
    apply_coverage_mask(var, dim_names, mask)

Set `var` to NaN wherever `mask` (over dimensions `dim_names`, as produced by
`coverage_mask`) says the observation has no data, broadcasting over time.

Why this is required: satellite retrievals do not cover the whole globe — MAC
`lwp` is ocean-only and is NaN over land, missing ~54% of grid points. Because
`zonal_average` ignores NaNs, the OBSERVED zonal mean is an average over covered
(ocean) points only, while the simulation has no NaNs and would be averaged over
ALL longitudes. That compares two different spatial samples, and for `lwp` it is
not a small effect: it flips the sign of the area-weighted bias (model 11.5% LOW
against the all-longitude sample vs 9.1% HIGH against the ocean-only sample,
because model LWP over land is much lower than over ocean — e.g. at 36.8°N the
all-longitude mean is 0.068 while the ocean-only mean is 0.116 against an observed
0.107). Masking the simulation to the observation's coverage first makes the two
zonal means like-for-like.

A no-op for variables whose observation covers every point (e.g. GPCP `pr`).
"""
function apply_coverage_mask(var, dim_names, mask)
    tname = ClimaAnalysis.time_name(var)
    var_names = filter(!=(tname), ordered_dim_names(var))
    # Compare on canonical names: obs and sim may spell the same axis differently
    # (CALIPSO `height` vs model `z`).
    var_canon = canonical_dim_name.(var_names)
    mask_canon = canonical_dim_name.(dim_names)
    if Set(var_canon) != Set(mask_canon)
        @warn "coverage mask dimensions do not match variable; skipping mask" short_name =
            ClimaAnalysis.short_name(var) var_names dim_names
        return var
    end
    # Reorder the mask axes to match this variable's axis order.
    perm = [findfirst(==(n), mask_canon) for n in var_canon]
    m = permutedims(mask, perm)

    data = copy(Array(var.data))
    if isnothing(tname) || !(tname in ordered_dim_names(var))
        size(m) == size(data) || (@warn "mask size mismatch; skipping"; return var)
        data[m] .= NaN
    else
        tidx = var.dim2index[tname]
        size(m) == size(selectdim(data, tidx, 1)) ||
            (@warn "mask size mismatch; skipping"; return var)
        for i in axes(data, tidx)
            selectdim(data, tidx, i)[m] .= NaN
        end
    end
    @info "Applied observational coverage mask to $(ClimaAnalysis.short_name(var))" masked_fraction =
        round(count(m) / length(m); digits = 3)
    return ClimaAnalysis.remake(var; data)
end

"""
    coarsen_lonlat(var, factor)

Reduce `var` to a coarse longitude-latitude grid by block averaging. Each
output cell is the mean of a `factor` x `factor` block of input cells, weighted
by `cosd(latitude)` within the block and skipping NaN cells. Block coordinates
are the unweighted centers of the input coordinates. Other dimensions (time,
altitude) pass through unchanged.

This is the aggregation counterpart of `zonal_average` for calibrations that
keep two spatial dimensions. `ClimaAnalysis.resampled_as` is not suitable here:
it interpolates point values at the coarse nodes, which discards most of the
data, aliases small scales, and spreads NaN from masked cells.

The longitude and latitude sizes must divide by `factor`. A NaN block (for
example a fully land-masked ocean retrieval block) stays NaN and is handled by
the flatten step like any other missing point.
"""
function coarsen_lonlat(var, factor)
    ClimaAnalysis.has_longitude(var) || return var
    lonname = ClimaAnalysis.longitude_name(var)
    latname = ClimaAnalysis.latitude_name(var)
    li = var.dim2index[lonname]
    lj = var.dim2index[latname]
    nlon = length(var.dims[lonname])
    nlat = length(var.dims[latname])
    (nlon % factor == 0 && nlat % factor == 0) || error(
        "coarsen_lonlat: grid $(nlon)x$(nlat) does not divide by factor $factor",
    )
    @info "Coarsening $(ClimaAnalysis.short_name(var)) by $factor: " *
          "$(nlon)x$(nlat) -> $(div(nlon, factor))x$(div(nlat, factor))"

    lats = collect(var.dims[latname])
    w_lat = cosd.(lats)

    outsize = collect(size(var.data))
    outsize[li] = div(nlon, factor)
    outsize[lj] = div(nlat, factor)
    out = fill(NaN, outsize...)

    for idx in CartesianIndices(Tuple(outsize))
        lonrange = ((idx[li] - 1) * factor + 1):(idx[li] * factor)
        latrange = ((idx[lj] - 1) * factor + 1):(idx[lj] * factor)
        num = 0.0
        den = 0.0
        for a in lonrange, b in latrange
            src = ntuple(
                d -> d == li ? a : d == lj ? b : idx[d], ndims(var.data),
            )
            x = var.data[src...]
            isfinite(x) || continue
            num += w_lat[b] * x
            den += w_lat[b]
        end
        den > 0 && (out[idx] = num / den)
    end

    blockmean(v) = [
        Statistics.mean(v[((k - 1) * factor + 1):(k * factor)]) for
        k in 1:div(length(v), factor)
    ]
    new_dims = copy(var.dims)
    new_dims[lonname] = blockmean(collect(var.dims[lonname]))
    new_dims[latname] = blockmean(lats)
    return ClimaAnalysis.remake(var; data = out, dims = new_dims)
end

"""
    reduce_spatial(var)

Apply the configured spatial reduction: `coarsen_lonlat` when the config
defines `COARSEN_FACTOR`, `zonal_average` otherwise. Must be applied
identically to the observations and the simulation so the flattened vectors
stay aligned.

A config may define `ZONAL_SHORT_NAMES` to reduce specific observables to
zonal means while the rest keep the coarsened 2-D grid. Use this to grade
an observable on its large-scale (amplitude) misfit when no calibrated
parameter can reach its spatial pattern; the full-field term otherwise
finances degradation of pattern-reachable observables (see the
lwp_clt_swcre_release verdict).
"""

"Collapse `(t0, t1)` sample/covariance ranges to `(t0, t0)`: after
`ClimaAnalysis.average_season_across_time` each season is one slice stamped
at its first date (SON -> Sep 1), so the builders select by that date."
collapse_ranges(ranges) = [(t0, t0) for (t0, _) in ranges]

function reduce_spatial(var)
    zonal_names =
        isdefined(Main, :ZONAL_SHORT_NAMES) ? Main.ZONAL_SHORT_NAMES : String[]
    if ClimaAnalysis.short_name(var) in zonal_names
        return zonal_average(var)
    end
    if isdefined(Main, :COARSEN_FACTOR)
        return coarsen_lonlat(var, Main.COARSEN_FACTOR)
    end
    return zonal_average(var)
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

    # Union of NaN locations across the selected time slices (spatial dims only),
    # shared with `coverage_mask` so the observation and the simulation are
    # guaranteed to be restricted to the same points.
    nanmask = date_range_nan_mask(var, date_ranges)

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
