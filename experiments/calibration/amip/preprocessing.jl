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

#Reduce `var` to a NaN-aware zonal (longitude) mean.
function zonal_average(var)
    ClimaAnalysis.has_longitude(var) || return var
    @info "Zonal (longitude) averaging $(ClimaAnalysis.short_name(var))"
    return ClimaAnalysis.average_lon(var; ignore_nan = true)
end

# Path of the saved observational coverage mask (see `coverage_mask` /
# `apply_coverage_mask`).
coverage_mask_path(output_dir) = joinpath(output_dir, "coverage_masks.jld2")

# Dimension names of `var` in the order of its data axes.
function ordered_dim_names(var)
    idx_name = sort([(idx, name) for (name, idx) in var.dim2index])
    return [name for (_, name) in idx_name]
end

"""
    date_range_nan_mask(var, date_ranges)

Union of `var`'s NaN points over the time slices in `date_ranges`: `true`
where the value is missing at any of those dates. Restricted to `date_ranges`
because a union over a decades-long record saturates to everything-missing.
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

Return `(dim_names, mask)` marking where `var` has no data over `date_ranges`.
`mask` is a `Bool` array over the non-time dimensions (`true` = missing) and
`dim_names` are those dimensions in mask axis order. `apply_coverage_mask`
uses it to sample the simulation at the same points as the observation.
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

Set `var` to NaN wherever `mask` (over `dim_names`, from `coverage_mask`) says
the observation has no data, broadcast over time.

Retrievals are not global (MAC `lwp` is ocean-only), so the observed zonal
mean covers only the observed points. Masking the simulation first makes the
two zonal means average over the same sample; without it the lwp bias flips
sign. A no-op when the observation covers every point.
"""
function apply_coverage_mask(var, dim_names, mask)
    tname = ClimaAnalysis.time_name(var)
    var_names = filter(!=(tname), ordered_dim_names(var))
    # Compare on conventional names: obs and sim may spell the same axis
    # differently (e.g. `height` vs `z`).
    var_canon = ClimaAnalysis.conventional_dim_name.(var_names)
    mask_canon = ClimaAnalysis.conventional_dim_name.(dim_names)
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
    block_mean(values, factor)

Mean of each consecutive block of `factor` entries of `values`. Trailing
entries that do not fill a block are dropped.
"""
block_mean(values, factor) = [
    Statistics.mean(values[((k - 1) * factor + 1):(k * factor)]) for
    k in 1:div(length(values), factor)
]

"""
    coarsen_lonlat(var, factor)

Block-average `var` onto a coarse longitude-latitude grid: each output cell
is the cosd(latitude)-weighted mean of a `factor` x `factor` block of input
cells, skipping NaNs. Block coordinates are the centers of the input
coordinates; other dimensions pass through. The lon and lat sizes must divide
by `factor`.

This is the 2D counterpart of `zonal_average` for calibrations that keep
geography. `ClimaAnalysis.resampled_as` is not suitable: it interpolates
point values at the coarse nodes, which discards most of the data and
spreads NaN from masked cells. A fully-NaN block stays NaN.
"""
function coarsen_lonlat(var, factor)
    ClimaAnalysis.has_longitude(var) || return var
    lonname = ClimaAnalysis.longitude_name(var)
    latname = ClimaAnalysis.latitude_name(var)
    # Extract longitude and latitude dimensions
    lon_axis = var.dim2index[lonname]
    lat_axis = var.dim2index[latname]
    nlon = length(var.dims[lonname])
    nlat = length(var.dims[latname])
    (nlon % factor == 0 && nlat % factor == 0) ||
        error("coarsen_lonlat: grid $(nlon)x$(nlat) does not divide by factor $factor")
    @info "Coarsening $(ClimaAnalysis.short_name(var)) by $factor: " *
          "$(nlon)x$(nlat) -> $(div(nlon, factor))x$(div(nlat, factor))"

    # Latitude cosine weighting
    lats = collect(var.dims[latname])
    lat_weights = cosd.(lats)

    # The output has the shape of the input with the two horizontal axes
    # divided by factor
    out_size = collect(size(var.data))
    out_size[lon_axis] = div(nlon, factor)
    out_size[lat_axis] = div(nlat, factor)
    out = fill(NaN, out_size...)

    for out_idx in CartesianIndices(Tuple(out_size))
        # Input cells that this output cell covers
        out_lon = out_idx[lon_axis]
        out_lat = out_idx[lat_axis]
        lon_block = ((out_lon - 1) * factor + 1):(out_lon * factor)
        lat_block = ((out_lat - 1) * factor + 1):(out_lat * factor)

        weighted_sum = 0.0
        weight_sum = 0.0
        for i_lon in lon_block, i_lat in lat_block
            # Follow out_idx on every axis except the two horizontal ones,
            # which walk the block
            src_idx = ntuple(ndims(var.data)) do d
                d == lon_axis ? i_lon : d == lat_axis ? i_lat : out_idx[d]
            end
            x = var.data[src_idx...]
            # Leave masked cells out of both sums
            isfinite(x) || continue
            weighted_sum += lat_weights[i_lat] * x
            weight_sum += lat_weights[i_lat]
        end
        # A fully masked block collects no weight and keeps its NaN
        weight_sum > 0 && (out[out_idx] = weighted_sum / weight_sum)
    end

    # A coarse coordinate is the mean of the input coordinates it covers
    new_dims = copy(var.dims)
    new_dims[lonname] = block_mean(collect(var.dims[lonname]), factor)
    new_dims[latname] = block_mean(lats, factor)
    return ClimaAnalysis.remake(var; data = out, dims = new_dims)
end

"""
    reduce_spatial(var, coarsen_factor)

Apply the configured spatial reduction: `coarsen_lonlat` by `coarsen_factor`,
or `zonal_average` when `coarsen_factor` is `nothing`. Apply identically to
obs and sim so the flattened vectors stay aligned.
"""
function reduce_spatial(var, coarsen_factor)
    isnothing(coarsen_factor) && return zonal_average(var)
    return coarsen_lonlat(var, coarsen_factor)
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

Force one shared NaN pattern across the dates in `date_ranges`: a point
missing at any of those dates becomes NaN at every time slice.

`SVDplusDCovariance` drops NaNs per sample and errors when the samples differ
in length; satellite coverage varies year to year, so each covariance date
would otherwise drop a different number of points.
"""
function harmonize_nan_mask_over_dates!(var, date_ranges)
    tname = ClimaAnalysis.time_name(var)
    isnothing(tname) && return var
    tidx = var.dim2index[tname]

    # Same mask construction as coverage_mask, so the observation and the
    # simulation are restricted to the same points.
    nanmask = date_range_nan_mask(var, date_ranges)

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
deviation for `var`. Skips NaNs.
"""
function compute_mean_and_stddev(var::ClimaAnalysis.OutputVar)
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
