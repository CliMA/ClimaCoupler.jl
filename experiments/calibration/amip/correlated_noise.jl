# Correlated diagonal noise for SVDplusD covariance
#
# The SVDplusD covariance is a low-rank interannual part plus a diagonal D. The
# diagonal treats every grid point as an independent constraint. On a 2-D grid
# points are strongly correlated, so a diagonal claims far more information
# than the field contains and the ensemble collapses.
#
# This file splits D into a downweighted diagonal and a spatially correlated
# term:
#
#     f * D + (1 - f) * sqrt(D) * K * sqrt(D),   K_ii = 1,   K_ij = exp(-d_ij / L)
#
# where d_ij is the great-circle distance between points i and j and L is the
# decorrelation length. f is the fraction of the floor variance that stays
# independent at each point, which is the retrieval noise and the
# representativeness error of a grid cell. The rest is model bias, which is
# spatially coherent. Points in different variables, levels, or times stay
# uncorrelated. Because K_ii = 1 the two terms sum to D on the diagonal, so the
# diagonal of the covariance is unchanged. A cluster of points inside one
# L-sized patch then counts as roughly one independent constraint.
#
# The result is stored as a dense `LinearAlgebra.Symmetric` matrix instead of an
# `EKP.SVDplusD` (whose floor must be `Diagonal`). Consumers that only read the
# covariance diagonal (steering indicators, preflight) are unaffected.

import LinearAlgebra
import ClimaAnalysis
import EnsembleKalmanProcesses as EKP

const EARTH_RADIUS = 6.371e6

"""
    great_circle_distance(lat1, lon1, cos_lat1, lat2, lon2, cos_lat2)

Great-circle distance in meters between two points, from latitudes and
longitudes in radians and the two latitude cosines. The caller precomputes the
trigonometry once per point instead of once per pair.
"""
function great_circle_distance(lat1, lon1, cos_lat1, lat2, lon2, cos_lat2)
    # This could use ClimaCore.Geometry.unit_great_circle_distance
    a = sin((lat2 - lat1) / 2)^2 + cos_lat1 * cos_lat2 * sin((lon2 - lon1) / 2)^2
    return 2 * EARTH_RADIUS * asin(min(1.0, sqrt(a)))
end

"""
    group_indices(group_keys)

Map each distinct entry of `group_keys` to the flattened indices that carry it,
so the pair loop can walk one level-time slice at a time.
"""
function group_indices(group_keys)
    groups = Dict{eltype(group_keys), Vector{Int}}()
    for (i, key) in enumerate(group_keys)
        push!(get!(() -> Int[], groups, key), i)
    end
    return groups
end

"""
    flattened_point_geometry(metadata)

Return `(lats, lons, group_keys)` for the points kept by the flattening
described by `metadata` (a `ClimaAnalysis` flatten `Metadata`), in flattened
order. `group_keys[i]` collects the non-horizontal coordinates (level, time) of
point `i`; only points with equal keys are correlated.
"""
function flattened_point_geometry(metadata)
    dims = collect(metadata.ordered_dims)
    conv = ClimaAnalysis.conventional_dim_name.(dims)
    # Coordinate values along each dimension, in flattening order
    axes_vals = [collect(metadata.dims[d]) for d in dims]
    sz = Tuple(length.(axes_vals))
    # Turns a linear index into the full grid into one index per dimension
    ci = CartesianIndices(sz)
    # Linear indices of the points that the flattening keeps
    kept = findall(!, metadata.drop_mask)

    lat_index = findfirst(==("latitude"), conv)
    isnothing(lat_index) && error("Correlated floor needs a latitude dimension; got $conv")
    lon_index = findfirst(==("longitude"), conv)

    # Positions of the non-horizontal dimensions (level and time)
    other_is = [i for i in eachindex(dims) if i != lat_index && i != lon_index]

    # For each kept point, read its latitude index out of ci, then its value
    lats = [axes_vals[lat_index][ci[k][lat_index]] for k in kept]
    # Assume zonal mean if there is no longitude dimension
    lons =
        isnothing(lon_index) ? zeros(length(kept)) :
        [axes_vals[lon_index][ci[k][lon_index]] for k in kept]
    # One tuple of level and time values per kept point
    group_keys = [Tuple(axes_vals[i][ci[k][i]] for i in other_is) for k in kept]
    return lats, lons, group_keys
end

"""
    add_correlated_floor!(
        total_cov, floor_sd, offset, metadata, decorrelation_length;
        identity_fraction,
    )

Add `f * D + (1 - f) * sqrt(D) * K * sqrt(D)` for one variable block into
`total_cov`, starting at row and column `offset + 1`, and return the number of
points added. `f` is `identity_fraction`, and `floor_sd` holds `sqrt(D)` for the
whole covariance.

`K` is the correlation kernel, with `d_ij` the great-circle distance in meters:

    K_ii = 1
    K_ij = exp(-d_ij / decorrelation_length)
    K_ij = 0     for i and j on different levels or at different times

`f` is the share of the floor variance that stays independent at each point.
Because `K_ii = 1`, the two terms sum to `D` on the diagonal, so `total_cov`
gains exactly `D` there. `f` also holds the smallest eigenvalue of the floor at
or above `f * D`, which keeps the covariance invertible at large
`decorrelation_length`, where the correlated term on its own approaches rank 1.
"""
function add_correlated_floor!(
    total_cov,
    floor_sd,
    offset,
    metadata,
    decorrelation_length;
    identity_fraction,
)
    lats, lons, group_keys = flattened_point_geometry(metadata)
    m = length(lats)
    # Per-point trigonometry, reused by every pair the point takes part in
    lats_rad = deg2rad.(lats)
    lons_rad = deg2rad.(lons)
    cos_lats = cos.(lats_rad)
    scale = 1.0 - identity_fraction

    # f * D from the first term plus (1 - f) * D from the second, since K_ii = 1
    for i in 1:m
        total_cov[offset + i, offset + i] += floor_sd[offset + i]^2
    end

    # Only points within one level/time slice correlate, so walk the pairs of
    # each slice
    for idxs in values(group_indices(group_keys))
        for (p, a) in enumerate(idxs), b in view(idxs, (p + 1):lastindex(idxs))
            d = great_circle_distance(
                lats_rad[a],
                lons_rad[a],
                cos_lats[a],
                lats_rad[b],
                lons_rad[b],
                cos_lats[b],
            )
            ga = offset + a
            gb = offset + b
            c = scale * exp(-d / decorrelation_length) * floor_sd[ga] * floor_sd[gb]
            total_cov[ga, gb] += c
            total_cov[gb, ga] += c
        end
    end
    return m
end

"""
    correlate_noise_floor(cov, metadata_vec, decorrelation_length; identity_fraction)

Split the diagonal floor D of an `EKP.SVDplusD` covariance into
`f * D + (1 - f) * sqrt(D) * K * sqrt(D)` and return the full covariance as a
dense `LinearAlgebra.Symmetric` matrix. `metadata_vec` holds one flatten
`Metadata` per variable, in the same order as the covariance blocks. The
covariance diagonal is unchanged.
"""
function correlate_noise_floor(
    cov::EKP.SVDplusD,
    metadata_vec,
    decorrelation_length;
    identity_fraction,
)
    n = first(size(cov))
    # Reconstruct the low-rank interannual part as a dense matrix, then
    # accumulate the floor blocks into it
    svd_cov = EKP.get_svd_cov(cov)
    total_cov = svd_cov.U * LinearAlgebra.Diagonal(svd_cov.S) * svd_cov.Vt
    floor_sd = sqrt.(EKP.get_diag_cov(cov).diag)

    off = 0
    for md in metadata_vec
        off += add_correlated_floor!(
            total_cov,
            floor_sd,
            off,
            md,
            decorrelation_length;
            identity_fraction,
        )
    end
    off == n || error("Metadata blocks cover $off points but the covariance has $n rows")
    return LinearAlgebra.Symmetric(total_cov)
end

"""
    apply_floor_correlation(obs_vec, decorrelation_length; identity_fraction)

Rebuild every `EKP.Observation` in `obs_vec` with a correlated noise floor.
All observations share one covariance.

This is the entry point, so it carries the only `identity_fraction` default.
The functions it calls require the value.
"""
function apply_floor_correlation(obs_vec, decorrelation_length; identity_fraction = 0.05)
    corr_cov = correlate_noise_floor(
        only(EKP.get_covs(first(obs_vec))),
        EKP.get_metadata(first(obs_vec)),
        decorrelation_length;
        identity_fraction,
    )
    return map(obs_vec) do obs
        EKP.Observation(
            Dict(
                "samples" => EKP.get_obs(obs),
                "covariances" => corr_cov,
                "names" => first(EKP.get_names(obs)),
                "metadata" => EKP.get_metadata(obs),
            ),
        )
    end
end
