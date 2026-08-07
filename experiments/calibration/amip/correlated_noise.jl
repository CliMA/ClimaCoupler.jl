# Correlated noise floor for the observation covariance.
#
# The SVDplusD covariance is a low-rank interannual part plus a diagonal floor
# D. The diagonal treats every grid point as an independent constraint. On a
# 2-D grid the points are strongly correlated, so a diagonal floor claims far
# more information than the field contains and the ensemble collapses (see the
# lwp_pr_2d10 and lwp_pr_2d15 runs).
#
# This file replaces the diagonal floor with a spatially correlated one:
#
#     floor = sqrt(D) * K * sqrt(D),   K_ij = exp(-d_ij / L)
#
# where d_ij is the great-circle distance between points i and j and L is the
# decorrelation length. Points in different variables, levels, or times stay
# uncorrelated. The diagonal of the covariance is unchanged because K has a
# unit diagonal. A cluster of points inside one L-sized patch then counts as
# roughly one independent constraint.
#
# The result is stored as a dense `LinearAlgebra.Symmetric` matrix instead of
# an `EKP.SVDplusD` (whose floor must be `Diagonal`). Consumers that only read
# the covariance diagonal (steering indicators, preflight) are unaffected.

import LinearAlgebra
import ClimaAnalysis
import EnsembleKalmanProcesses as EKP

const EARTH_RADIUS_M = 6.371e6

"""
    great_circle_distance_m(lat1, lon1, lat2, lon2)

Great-circle distance in meters between two points given in degrees.
"""
function great_circle_distance_m(lat1, lon1, lat2, lon2)
    phi1 = deg2rad(lat1)
    phi2 = deg2rad(lat2)
    dphi = phi2 - phi1
    dlam = deg2rad(lon2 - lon1)
    a = sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlam / 2)^2
    return 2 * EARTH_RADIUS_M * asin(min(1.0, sqrt(a)))
end

"""
    flattened_point_geometry(metadata)

Return `(lats, lons, group_keys)` for the points kept by the flattening
described by `metadata` (a `ClimaAnalysis` flatten `Metadata`), in flattened
order. `lons` is `nothing` when the variable has no longitude dimension (for
example after a zonal mean). `group_keys[i]` collects the non-horizontal
coordinates (level, time) of point `i`; only points with equal keys are
correlated.
"""
function flattened_point_geometry(metadata)
    dims = collect(metadata.ordered_dims)
    conv = ClimaAnalysis.conventional_dim_name.(dims)
    axes_vals = [collect(metadata.dims[d]) for d in dims]
    sz = Tuple(length.(axes_vals))
    ci = CartesianIndices(sz)
    kept = findall(!, metadata.drop_mask)

    lat_i = findfirst(==("latitude"), conv)
    isnothing(lat_i) && error("Correlated floor needs a latitude dimension; got $conv")
    lon_i = findfirst(==("longitude"), conv)
    other_is = [i for i in eachindex(dims) if i != lat_i && i != lon_i]

    lats = [axes_vals[lat_i][ci[k][lat_i]] for k in kept]
    lons = isnothing(lon_i) ? nothing : [axes_vals[lon_i][ci[k][lon_i]] for k in kept]
    group_keys = [Tuple(axes_vals[i][ci[k][i]] for i in other_is) for k in kept]
    return lats, lons, group_keys
end

"""
    floor_correlation_block(metadata, decorrelation_length; identity_fraction)

Dense correlation matrix for one variable block. Entry `(i, j)` is
`exp(-d_ij / decorrelation_length)` when points `i` and `j` share the same
level and time, and `0` otherwise. The result is blended with the identity
(`identity_fraction`) to guarantee a well-conditioned kernel; the diagonal
stays exactly 1.
"""
function floor_correlation_block(metadata, decorrelation_length; identity_fraction = 0.05)
    lats, lons, group_keys = flattened_point_geometry(metadata)
    n = length(lats)
    K = Matrix{Float64}(LinearAlgebra.I, n, n)
    scale = 1.0 - identity_fraction
    for a in 1:n, b in (a + 1):n
        group_keys[a] == group_keys[b] || continue
        d = if isnothing(lons)
            EARTH_RADIUS_M * deg2rad(abs(lats[a] - lats[b]))
        else
            great_circle_distance_m(lats[a], lons[a], lats[b], lons[b])
        end
        c = scale * exp(-d / decorrelation_length)
        K[a, b] = c
        K[b, a] = c
    end
    return K
end

"""
    correlate_noise_floor(cov, metadata_vec, decorrelation_length; identity_fraction)

Turn the diagonal floor D of an `EKP.SVDplusD` covariance into a spatially
correlated floor `sqrt(D) * K * sqrt(D)` and return the full covariance as a
dense `LinearAlgebra.Symmetric` matrix. `metadata_vec` holds one flatten
`Metadata` per variable, in the same order as the covariance blocks. The
covariance diagonal is unchanged.
"""
function correlate_noise_floor(
    cov::EKP.SVDplusD,
    metadata_vec,
    decorrelation_length;
    identity_fraction = 0.05,
)
    n = first(size(cov))
    svdc = EKP.get_svd_cov(cov)
    sigma = svdc.U * LinearAlgebra.Diagonal(svdc.S) * svdc.Vt
    sd = sqrt.(EKP.get_diag_cov(cov).diag)

    off = 0
    for md in metadata_vec
        K = floor_correlation_block(md, decorrelation_length; identity_fraction)
        m = size(K, 1)
        rng = (off + 1):(off + m)
        sigma[rng, rng] .+= (sd[rng] * sd[rng]') .* K
        off += m
    end
    off == n || error("Metadata blocks cover $off points but the covariance has $n rows")
    return LinearAlgebra.Symmetric(sigma)
end

"""
    apply_floor_correlation(obs_vec, decorrelation_length; identity_fraction)

Rebuild every `EKP.Observation` in `obs_vec` with the correlated noise floor.
All observations share one covariance (it is estimated from the covariance
dates, not from the sample), so the dense matrix is computed once and reused.
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
