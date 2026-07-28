# Steering indicators for a running calibration.
#
# `steering_summary` builds a plain-language health block for one iteration.
# The caller logs the block. Nothing is written to disk. The indicators are
# advisory. They do not stop or change the run.
#
# The six indicators:
#   1. Fit-to-noise for each observable.
#   2. Ensemble collapse and contraction.
#   3. Parameter drift and a convergence hint.
#   4. Identifiability of each parameter.
#   5. Reachability of each observable.
#   6. Degeneracy between parameter pairs.
#
# All quantities come from the ekp object, the current G ensemble, and the
# prior. The ekp object holds the full parameter history, so trend indicators
# need no external state.

import Statistics as SIStatistics
import EnsembleKalmanProcesses as SIEKP

# All thresholds in one table.
const STEERING_THRESHOLDS = (;
    # RMS residual per observable, in units of the observation noise.
    obs_signal = 2.0,          # above: learnable signal remains
    obs_floor_lo = 0.7,        # below: overfit risk, or the noise is too large
    # Ensemble spread for each parameter, as a ratio to the prior spread.
    collapse = 1e-3,           # below: the ensemble has collapsed
    no_learning = 0.9,         # near 1 after several iterations: no learning
    # Mean parameter movement per iteration, in units of the current ensemble
    # spread. A converged ensemble still wanders inside its posterior, so prior
    # units would keep the movement above any small threshold forever.
    converged_drift = 0.15,
    # Fraction of observation points inside the ensemble G envelope.
    envelope_ok = 0.5,
    onesided = 0.9,            # sign consistency of the residual
    # Correlation limit for a parameter pair.
    degenerate = 0.9,
)

"""
    _obs_blocks(ekp, iteration)

Return `(y, sigma, blocks)` for the current minibatch. `blocks` holds the name
and index range of each observed variable in the flattened observation vector.
"""
function _obs_blocks(ekp, iteration)
    obs_series = SIEKP.get_observation_series(ekp)
    minibatch_obs =
        ClimaCalibrate.get_observations_for_nth_iteration(obs_series, iteration)
    y = mapreduce(SIEKP.get_obs, vcat, minibatch_obs)
    sigma = mapreduce(vcat, minibatch_obs) do obs
        cov = SIEKP.get_obs_noise_cov(obs)
        [sqrt(abs(cov[i, i])) for i in 1:size(cov, 1)]
    end
    obs_vars = mapreduce(
        ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, minibatch_obs,
    )
    names = ClimaAnalysis.short_name.(obs_vars)
    lens = [length(ClimaAnalysis.flatten(v).data) for v in obs_vars]
    blocks = []
    off = 0
    for (n, l) in zip(names, lens)
        push!(blocks, (; name = n, rng = (off + 1):min(off + l, length(y))))
        off += l
    end
    return y, sigma, blocks
end

"Row means of the G ensemble. NaN entries are skipped."
function _g_mean(g_ensemble)
    return map(eachrow(g_ensemble)) do row
        finite = filter(isfinite, row)
        isempty(finite) ? NaN : SIStatistics.mean(finite)
    end
end

"""
    steering_summary(ekp, g_ensemble, prior, iteration) -> String

Build the steering block for this iteration. The function is pure. Wrap the
call in try/catch at the call site.
"""
function steering_summary(ekp, g_ensemble, prior, iteration)
    T = STEERING_THRESHOLDS
    lines = String["── steering (iter $iteration) " * "─"^40]

    # Indicators 1 and 5: fit-to-noise and reachability for each observable.
    y, sigma, blocks = _obs_blocks(ekp, iteration)
    gm = _g_mean(g_ensemble)
    n = min(length(gm), length(y), length(sigma))
    obs_done = Bool[]
    for b in blocks
        idx = [i for i in b.rng if i <= n && isfinite(gm[i]) && sigma[i] > 0]
        isempty(idx) && continue
        resid = [(gm[i] - y[i]) / sigma[i] for i in idx]
        rms = sqrt(SIStatistics.mean(abs2, resid))
        verdict = rms > T.obs_signal ? "learnable signal remains" :
                  rms < T.obs_floor_lo ?
                  "below noise floor. Overfit risk, or the noise is too large" :
                  "fit to noise floor"
        push!(obs_done, T.obs_floor_lo <= rms <= 1.3)
        push!(lines, "obs    " * rpad(b.name, 6) * " RMS " *
                     rpad(string(round(rms; digits = 2)) * "σ", 7) * " ▸ " * verdict)

        inenv = 0
        tot = 0
        sgn = 0.0
        for i in idx
            col = view(g_ensemble, i, :)
            finite = filter(isfinite, col)
            isempty(finite) && continue
            tot += 1
            minimum(finite) <= y[i] <= maximum(finite) && (inenv += 1)
            sgn += sign(gm[i] - y[i])
        end
        if tot > 0
            frac = inenv / tot
            onesided = abs(sgn / tot)
            flag = (frac < T.envelope_ok && onesided > T.onesided && rms > T.obs_signal) ?
                   " ▸ STRUCTURAL. No parameter direction reaches this observable. Raise its noise floor or remove it" :
                   frac >= T.envelope_ok ? " ▸ reachable" : " ▸ partially reachable"
            push!(lines, "reach  " * rpad(b.name, 6) * " " *
                         string(round(Int, 100 * frac)) * "% in envelope" * flag)
        end
    end

    # Indicators 2 and 4: collapse, contraction, and identifiability.
    pnames = SIEKP.get_name(prior)
    u_hist = SIEKP.get_u(ekp)
    u0, unow = first(u_hist), last(u_hist)
    spread0 = vec(SIStatistics.std(u0; dims = 2))
    spreadn = vec(SIStatistics.std(unow; dims = 2))
    ratios = [s0 > 0 ? sn / s0 : NaN for (sn, s0) in zip(spreadn, spread0)]
    collapsed = any(r -> isfinite(r) && r < T.collapse, ratios)
    stuck = all(r -> isfinite(r) && r > T.no_learning, ratios) && length(u_hist) > 3
    spread_str = join(
        ["$(pn) $(round(r; digits = 2))×prior" for (pn, r) in zip(pnames, ratios)],
        " · ",
    )
    spread_verdict = collapsed ?
        "COLLAPSED. The observation covariance is over-informative. Stop and fix the noise" :
        stuck ? "no contraction. The parameter signal is below the noise" :
        "contracting, no collapse"
    push!(lines, "spread " * spread_str * " ▸ " * spread_verdict)
    ident = join(
        [r > T.no_learning && length(u_hist) > 3 ? "$(pn) unconstrained" :
         "$(pn) constrained" for (pn, r) in zip(pnames, ratios)],
        " · ",
    )
    push!(lines, "ident  " * ident *
                 (any(r -> r > T.no_learning, ratios) && length(u_hist) > 3 ?
                  " ▸ unconstrained parameters carry no data signal. Add an observable or remove them" : ""))

    # Indicator 6: degeneracy between parameter pairs.
    phinow = SIEKP.get_ϕ(prior, ekp, length(u_hist))
    p = size(phinow, 1)
    degen_msgs = String[]
    for i in 1:p, j in (i + 1):p
        c = SIStatistics.cor(phinow[i, :], phinow[j, :])
        if isfinite(c) && abs(c) > T.degenerate
            push!(degen_msgs,
                  "|r($(pnames[i]), $(pnames[j]))| = $(round(abs(c); digits = 2)) ▸ NOT separately identifiable")
        end
    end
    push!(lines, isempty(degen_msgs) ? "degen  all parameter pairs separable" :
                 "degen  " * join(degen_msgs, " · "))

    # Indicator 3: parameter drift and the convergence hint. Drift is measured
    # in units of the current ensemble spread, so movement inside the settled
    # posterior does not block the convergence call.
    if length(u_hist) >= 2
        phiprev = SIEKP.get_ϕ(prior, ekp, length(u_hist) - 1)
        now_sigma = vec(SIStatistics.std(phinow; dims = 2))
        drift = [
            now_sigma[i] > 0 ?
            abs(SIStatistics.mean(phinow[i, :]) - SIStatistics.mean(phiprev[i, :])) /
            now_sigma[i] : NaN for i in 1:p
        ]
        drift_str = join(
            ["$(pn) $(round(d; digits = 2))×spread/iter" for (pn, d) in zip(pnames, drift)],
            " · ",
        )
        converged = all(d -> isfinite(d) && d < T.converged_drift, drift) &&
                    !isempty(obs_done) && all(obs_done)
        push!(lines, "drift  " * drift_str *
                     (converged ?
                      " ▸ CONVERGED candidate. More iterations are unlikely to learn more" :
                      ""))
    end

    push!(lines, "─"^58)
    return join(lines, "\n")
end
