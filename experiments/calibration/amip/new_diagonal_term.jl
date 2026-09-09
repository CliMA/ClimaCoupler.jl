import ClimaCalibrate
import ClimaCalibrate.ObservationRecipe as ObservationRecipe
import ClimaCalibrate.ObservationRecipe:
    AbstractDiagonalTerm, VarianceDiagonal, compute_diagonal
import EnsembleKalmanProcesses as EKP
import LinearAlgebra: Diagonal
import Statistics

"""
    ImprovedDiagonal <: AbstractDiagonalTerm

Diagonal term in the form of D_y + β * D_{diff} + σ^2 * s^2 * I where β, σ, and
s^2 are per variable and s^2 is the median of the variable's block of D_y (the
sample variances).
"""
struct ImprovedDiagonal <: AbstractDiagonalTerm
    betas::Vector{Float64}
    sigma_squares::Vector{Float64}
    rank::Union{Nothing, Int}
end

"""
    ImprovedDiagonal(betas = 0.05, sigma_squares = 10^-6; rank = nothing)

Constructor for the `ImprovedDiagonal`.

The `i`th entry of `betas` and `sigma_squares` is used for the `i`th variable.
If there is only one entry, then the same value is used for all variables.
"""
function ImprovedDiagonal(betas = 0.05, sigma_squares = 10^-6; rank = nothing)
    all(>=(0), sigma_squares) ||
        error("Sigma squares ($sigma_squares) should not be negative")
    return ImprovedDiagonal(
        float.(vcat(betas)),
        float.(vcat(sigma_squares)),
        rank,
    )
end

"""
    compute_diagonal(term::ImprovedDiagonal, sample_collection)

Compute the diagonal from `term` and `sample_collection`, a matrix of samples.
"""
function compute_diagonal(term::ImprovedDiagonal, sample_collection)
    S = ClimaCalibrate.SampleBuilder.get_samples(sample_collection)
    ranges = ClimaCalibrate.SampleBuilder.var_indices(sample_collection)
    # Check the number of β and σ^2 are what we expect
    for p in (term.betas, term.sigma_squares)
        length(p) in (1, length(ranges)) || error(
            "Expected 1 or $(length(ranges)) parameters per variable, got $p",
        )
    end

    # Compute D_y
    V_y = compute_diagonal(VarianceDiagonal(), sample_collection).diag
    gamma = isnothing(term.rank) ? EKP.tsvd_cov_from_samples(S) :
            EKP.tsvd_cov_from_samples(S, term.rank)
    gamma_diag = vec((gamma.U .^ 2) * gamma.S)
    diagonal = max.(V_y .- gamma_diag, 0)

    for (v, rng) in enumerate(ranges)
        block = view(S, rng, :)
        # s^2 is the median of V_y (the sample variances) for this variable
        s_square = Statistics.median(view(V_y, rng))
        # Add β * D_{diff} and σ^2 * s^2 * I
        # Q: Is it beta^2 or beta? What is a good value for beta?
        # D_{diff} is the per-entry mean over samples of the squared values
        # (row-wise mean of the block, one value per entry of the variable)
        d_diff = vec(Statistics.mean(abs2, block; dims = 2))
        diagonal[rng] .+=
            term.betas[v] .* d_diff .+
            term.sigma_squares[v] * s_square
    end
    all(>(0), diagonal) || error("Negative diagonal entries; check the sign of betas")
    return Diagonal(diagonal)
end

"""
    ImprovedSVDplusDCovariance{S <: ObservationRecipe.SVDplusDCovariance}

An observation recipe for a `SVDplusD` covariance matrix with a
`ImprovedDiagonal` matrix.
"""
struct ImprovedSVDplusDCovariance{S <: ObservationRecipe.SVDplusDCovariance} <:
       ObservationRecipe.AbstractCovarianceEstimator
    inner::S
end

"""
    ImprovedSVDplusDCovariance(;
        betas = 0.05,
        sigma_squares = 10^-6;
        rank = nothing,
        use_latitude_weights = false,
        min_cosd_lat = 0.1,
    )

Create a `SVDplusD` covariance matrix with a `ImprovedDiagonal` matrix.
"""
function ImprovedSVDplusDCovariance(;
    betas = 0.05^2, # TODO: Find a good value for beta
    sigma_squares = 10^-6;
    rank = 2,
    use_latitude_weights = true,
    min_cosd_lat = 0.1,
)
    return ImprovedSVDplusDCovariance(
        ObservationRecipe.SVDplusDCovariance(
            ImprovedDiagonal(betas, sigma_squares; rank);
            rank,
            use_latitude_weights,
            use_weight_samples_for_diagonal = use_latitude_weights,
            min_cosd_lat,
        ),
    )
end

"""
    ObservationRecipe.covariance(est::ImprovedSVDplusDCovariance, sample_collection)

Compute the covariance matrix from `est` and `sample_collection`, a matrix of
samples.
"""
ObservationRecipe.covariance(est::ImprovedSVDplusDCovariance, sample_collection) =
    ObservationRecipe.covariance(est.inner, sample_collection)
