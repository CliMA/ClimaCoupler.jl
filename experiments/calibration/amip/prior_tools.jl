# Shared prior constructors for calibration configs.

import EnsembleKalmanProcesses.ParameterDistributions as PriorToolsPD
import Random as PriorToolsRandom
import Statistics as PriorToolsStatistics

"""
    checked_constrained_gaussian(name, mu, sigma, lo, hi)

`constrained_gaussian` that stays correct for small-magnitude parameters.

EKP's `constrained_gaussian` fits the unconstrained normal by Nelder-Mead
from a unit-normal start, and the optimizer declares convergence on an
ABSOLUTE objective tolerance. For targets with mu, sigma below ~1e-2 the
objective is already under the tolerance at the start point, so the solver
returns the unit normal unchanged: the prior silently becomes a logistic
unit normal over the bounds (mean at the bounds center) and the requested
moments are ignored. The post-fit check also uses absolute tolerances, so
no warning is raised. This silently shaped every q_liq prior of the July
2026 campaign (mean 5e-4 or 7.5e-4 = bounds center, not the configured
value).

Fix: solve in sigma-scaled units, where the objective is O(1), and rebuild
on the original bounds. The logistic transform maps the same unconstrained
normal to both problems, so the solution transfers exactly. The result is
verified by sampling; construction fails loudly on a moment mismatch.
"""
function checked_constrained_gaussian(name, mu, sigma, lo, hi)
    (isfinite(lo) && isfinite(hi)) ||
        error("checked_constrained_gaussian requires finite bounds")
    scaled =
        PriorToolsPD.constrained_gaussian(name, mu / sigma, 1.0, lo / sigma, hi / sigma)
    normal = scaled.distribution[1].distribution
    d = PriorToolsPD.ParameterDistribution(
        PriorToolsPD.Parameterized(normal),
        PriorToolsPD.bounded(lo, hi),
        name,
    )
    c = PriorToolsPD.transform_unconstrained_to_constrained(
        d,
        PriorToolsPD.sample(PriorToolsRandom.MersenneTwister(0), d, 100_000),
    )
    m = PriorToolsStatistics.mean(c)
    s = PriorToolsStatistics.std(c)
    # atol scaled to sigma: with a zero-mean target (e.g. the Pi-groups c2/c3
    # priors) a pure rtol check can never pass - the tolerance collapses to
    # 0.05*|m| - even though a mean within a few percent of a sigma is an
    # excellent fit. For |mu| >> sigma the rtol term still dominates. The
    # silent-unit-normal failure this function exists to catch puts the mean
    # at the bounds center, which IS mu for symmetric zero-mean priors, so
    # there the std check below is the real guard either way.
    isapprox(m, mu; rtol = 0.05, atol = 0.05 * sigma) ||
        error("prior for $name missed its mean: target $mu, got $m")
    isapprox(s, sigma; rtol = 0.15) ||
        error("prior for $name missed its std: target $sigma, got $s")
    return d
end
