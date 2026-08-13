# Steering Indicators: per-iteration calibration health, logged in plain language

## Context

Every finding that steered this study — collapse, noise-dominated cl, "fit to within
noise," loose rain_τ, unreachable swcre — was computed manually, after the fact.
This adds a **steering block to driver.log every iteration** so the calibration
self-reports its health as it runs. **Log-only**: no files are written; nothing is
persisted beyond what ClimaCalibrate already saves. Advisory, never auto-acting.

All indicators are computable in-memory from what `analyze_iteration` already
receives — and the `ekp` object carries its full history (`get_ϕ` across iterations,
error series), so trend indicators need no external state.

## v1 — the first PR

**Six indicators, one log block, no files, no new dependencies, no figures.**

### Indicators (each with the incident that justifies it)

1. **Per-observable fit-to-noise** — RMS of `(mean(G) − y)/σ` per variable block
   (blocks via `ObservationRecipe.reconstruct_vars`; the number is already computed
   for the g_vs_obs plot — this records the verdict in the log).
   `>2σ` "learnable signal remains" · `≈1σ` "fit to noise floor — no further
   information from this observable" · `<0.7σ` "below noise floor — overfit risk or
   σ overestimated."
   *lwp went 2.22σ → 1.0σ in 3 iterations; that was the convergence signal.*

2. **Collapse / contraction monitor** — per-parameter ensemble spread
   (unconstrained space) as a ratio to prior spread.
   `<1e-3` "COLLAPSED — observation covariance over-informative" · stuck at ~1
   "no learning — parameter signal below noise."
   *The lwp+cl collapse: spread → 1e-13, parameters frozen to 16 digits.*

3. **Convergence hint** — per-parameter movement per iteration in prior-σ units
   (from the ekp's own ϕ history). All params `<0.05σ/iter` AND all obs ≤~1σ →
   "CONVERGED candidate — further iterations likely spend GPU without learning."
   *LWP-only ran 9 iterations; parameters settled by ~5.*

4. **Per-parameter identifiability** — posterior/prior variance ratio per
   parameter. Still ≈1 after several iterations → "data does not constrain this
   parameter with the current observables."
   *rain_τ stayed loose (1450–1670) under lwp-only; only pr pinned it.*

5. **Reachability / structural-bias flag** — per obs block: fraction of points
   where y lies inside the ensemble-G envelope, plus sign-consistency of the
   residual. Persistently one-sided >2σ block → "STRUCTURAL — no parameter
   direction reaches this observable; raise its noise floor or drop it, else it
   will distort reachable parameters."
   *swcre: obs −40 W/m² vs model (−69, −53) across the entire swept range.*

6. **Degeneracy monitor** — parameter-parameter correlation from the constrained
   ensemble (a p×p matrix, trivial at these sizes). |r|>0.9 → "these parameters are
   not separately identifiable with the current observables."
   *q_liq/rain_τ against lwp alone.*

### The log block (the entire output surface)

`@info`'d into driver.log each iteration from `analyze_iteration`:

```
── steering (iter 4) ────────────────────────────────────────
obs    lwp   RMS 1.02σ  ▸ fit to noise floor
obs    pr    RMS 0.98σ  ▸ fit to noise floor
reach  lwp   96% in envelope · pr 93% ▸ reachable
spread q_liq 0.31×prior · rain_τ 0.55×prior ▸ contracting, no collapse
ident  q_liq var 0.09×prior ▸ constrained · rain_τ 0.30×prior ▸ constrained
degen  |r(q_liq, rain_τ)| = 0.41 ▸ separable
drift  q_liq 0.03σ/iter · rain_τ 0.04σ/iter ▸ CONVERGED candidate
─────────────────────────────────────────────────────────────
```

### Files

- **New:** `experiments/calibration/amip/steering_indicators.jl` —
  `steering_summary(ekp, g_ensemble, prior, iteration) -> String` (pure; builds the
  block), thresholds in one `const STEERING_THRESHOLDS` table with a one-line "why"
  per entry. Experiments-local (promotion to src/CalibrationTools is a later
  refactor, keeping this PR free of package/test surface).
- **Modified:** `experiments/calibration/amip/post_analyze_iteration.jl` — ~5 lines
  in `analyze_iteration`: include + `@info steering_summary(...)`, wrapped in
  try/catch like the existing plot hooks (never fatal).

### Validation (acceptance test for the PR)

Retrospective script `test_indicators_retrospective.jl`: load saved eki/G files from
two existing run directories, print the steering block per iteration, and assert the
verdict strings:
1. `amip_calibration_out_collapsed_v1` → COLLAPSED flagged at iteration 2.
2. `amip_calibration_lwp_pr_out` → lwp 2.2σ→1.0σ; both obs reach "fit to noise
   floor"; CONVERGED candidate by ~iteration 4; no collapse; parameters separable.
Both narratives were verified by hand this study; the indicators must re-discover
them from the saved state alone.

## Deferred (backlog)

- Observation-influence audit (per-block loss contribution → noise-dominated /
  over-informing flags — the cl and collapse lessons at finer grain).
- Physical guards (global-mean pr, TOA imbalance, negative lwp per member).
- Operational health (member failures, NaN fraction, wall-time per iteration).
- Dashboard figure, report/sweep integration, STOP_RECOMMENDED marker,
  promotion to `src/CalibrationTools/Indicators.jl` with package tests.
