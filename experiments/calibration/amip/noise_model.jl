# The observation noise model used by the AMIP calibrations:
#
#     Gamma = C_y + D_y + beta * D_mean + sigma2 * s^2 * I
#
#   C_y      rank-`rank` truncated-SVD sample covariance of the observation years, supplied
#            by `SVDplusDCovariance`. Only one or two modes are reliably estimated from
#            O(10) samples, so the rest is not kept as low-rank structure but folded into
#   D_y      the per-point sample variance MINUS the variance the retained modes already
#            explain. Without this subtraction the diagonal of `C_y + D` counts the
#            retained modes twice. This is the only piece ClimaCalibrate does not provide.
#   D_mean   the square of the per-point sample mean, from `ModelErrorScaleDiagonal`.
#            `beta` sets how much of the model-data difference is treated as irreducible
#            model error, and is a VARIANCE scale: beta = 0.05^2 is a 5% model error.
#   s^2      the per-variable median of the positive sample variances, so the absolute
#            floor sigma2 * s^2 is scale aware. This term is there to keep the covariance
#            well conditioned, not to express a belief about uncertainty.
#
# Latitude weighting is applied to the samples before both the SVD and the diagonal, so
# that C_y and D_y are computed from the same matrix and the subtraction in D_y is valid.
#
# Everything here is a lazy description. ClimaCalibrate computes the matrix from the
# (weighted) SampleCollection when the observation is built.

import Statistics
import LinearAlgebra: Diagonal, svd
import ClimaCalibrate: SampleBuilder
import ClimaCalibrate.ObservationRecipe as ObservationRecipe

"""
    ResidualVarianceDiagonal(rank)

The per-point sample variance that the leading `rank` modes of the sample covariance do
not explain.

Use the same `rank` as the accompanying `SVDplusDCovariance`, so that the low-rank part
plus this diagonal reproduces the full sample variance. With `rank = 0` this is the full
variance, the same as `ObservationRecipe.VarianceDiagonal()`.
"""
struct ResidualVarianceDiagonal <: ObservationRecipe.AbstractDiagonalTerm
    rank::Int
    function ResidualVarianceDiagonal(rank::Integer)
        rank >= 0 || throw(ArgumentError("rank ($rank) must be non-negative"))
        return new(rank)
    end
end

function ObservationRecipe.compute_diagonal(
    term::ResidualVarianceDiagonal,
    sample_collection,
)
    samples = SampleBuilder.get_samples(sample_collection)
    FT = eltype(samples)
    n_samples = size(samples, 2)
    n_samples >= 2 || throw(
        ArgumentError(
            "ResidualVarianceDiagonal needs at least 2 samples, got $n_samples. " *
            "Give the observation more than one date range.",
        ),
    )
    # Columns of a scaled, centred sample matrix, so that sum(abs2) over columns is the
    # 1/(N-1) sample variance and the SVD of this matrix has the covariance eigenvalues
    # as its squared singular values.
    centred = (samples .- Statistics.mean(samples, dims = 2)) ./ sqrt(FT(n_samples - 1))
    total = vec(sum(abs2, centred, dims = 2))
    rank = min(term.rank, n_samples - 1)
    explained = if rank == 0
        zero(total)
    else
        factorization = svd(centred)
        vec(sum(abs2, factorization.U[:, 1:rank] .* factorization.S[1:rank]', dims = 2))
    end
    # max() guards round-off only; the subtraction is exact in theory.
    return Diagonal(convert(Vector{FT}, max.(total .- explained, zero(FT))))
end

"""
    noise_diagonal(; beta, rank, sigma2 = 1e-6)

The diagonal of `Gamma`: `D_y + beta * D_mean + sigma2 * s^2 * I`.

`beta` is a variance scale and may be a scalar or one value per variable, in the order the
variables appear in the `SampleCollection`. Pair this with an `SVDplusDCovariance` of the
same `rank`.
"""
function noise_diagonal(; beta, rank, sigma2 = 1e-6)
    betas = beta isa AbstractVector ? Float64.(collect(beta)) : [Float64(beta)]
    all(>=(0), betas) || throw(ArgumentError("beta must be non-negative, got $betas"))
    sigma2 >= 0 || throw(ArgumentError("sigma2 must be non-negative, got $sigma2"))
    return ResidualVarianceDiagonal(rank) .+
           ObservationRecipe.ModelErrorScaleDiagonal(sqrt.(betas)) .+
           sigma2 .*
           ObservationRecipe.QuantileDiagonal(0.5, ObservationRecipe.VarianceDiagonal())
end

"""
    noise_covariance_estimator(; beta, rank = 2, sigma2 = 1e-6,
                               use_latitude_weights = true, min_cosd_lat = 0.1)

An `SVDplusDCovariance` implementing `Gamma`, with latitude weighting applied to the
samples before both the SVD and the diagonal terms.
"""
function noise_covariance_estimator(;
    beta,
    rank = 2,
    sigma2 = 1e-6,
    use_latitude_weights = true,
    min_cosd_lat = 0.1,
)
    return ObservationRecipe.SVDplusDCovariance(
        noise_diagonal(; beta, rank, sigma2);
        rank,
        use_latitude_weights,
        use_weighted_samples_for_diagonal = true,
        min_cosd_lat,
    )
end
