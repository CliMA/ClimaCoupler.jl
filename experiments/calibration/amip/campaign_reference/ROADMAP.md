# Calibration Roadmap: from 2-parameter runs to the large joint calibration

## Context and objective

The end goal is **large-scale calibration: many parameters × many observations at
once**. The work so far (2-parameter warm-rain calibrations) is the derisking path to
that goal, not an alternative to it. Every failure mode a large calibration has was
already hit — cheaply — on the small ones:

| failure mode | where we hit it | defense now in hand |
|---|---|---|
| ensemble collapse from over-informative obs | lwp+cl, 35k constraints | zonal aggregation to true DOF, SVDplusD, model-error floor 0.2 |
| compensating errors from unreachable obs | swcre scouting (r=0.95 with q_liq, unreachable) | reachability gate + per-obs σ-residuals |
| sampling inconsistency | ocean-only lwp vs all-longitude sim (bias sign flip) | coverage masking (committed 6ffd3646) |
| silent environment drift | Atmos 0.41.3→0.42.2 rebase (27% LWP shift) | packaging/freeze framework (PACKAGING_PLAN.md) |
| overfitting one sample | single-Oct targets | multi-year minibatch machinery (SON runs) |

Central insight from the swcre/cl scouting: **degeneracy and reachability are
properties of the parameter set, not the observation set.** swcre is "a duplicate of
lwp" only in a 2-parameter space; with cloud-fraction and brightness parameters in
play it becomes the constraint that pins them. More parameters need more observations
to stay identifiable; more observations need more parameters to be reachable. They
must scale together — via the identifiability map.

## Current state (2026-07-28)

- **Validated template**: lwp+pr joint calibration on corrected setup (Atmos 0.42.2,
  coverage-masked obs), q_liq → ~3.1e-4, rain_τ → ~1150-1370, both obs fit to ~1σ,
  healthy spread. (Old q_liq≈7e-4 obsolete: masking bug + model upgrade.)
- **Measured noise floors** (zonal-mean, Oct, interannual + 20% floor):
  lwp 0.026 kg/m², pr 0.73 mm/day, swcre 9.8 W/m², lwcre 5.1 W/m², cl 0.043.
- **Known structural biases at best-fit warm-rain params**: cl ~18% too low ("too
  few"), swcre 13-28 W/m² too reflective ("too bright"), lwcre ~4% (nearly unbiased).
- **Warm-rain parameter family exhausted**: lwp/pr at noise floor; swcre/cl
  unreachable by q_liq/rain_τ.

## Phase 0 — close out current runs (days)

- Let lwp+pr finish; generate its calibration report; record posterior (q_liq,
  rain_τ) as the warm-rain result for Atmos 0.42.2.
- Let the SON spin-up-vs-no-spin-up diagnosis finish (internally valid differential
  comparison despite pre-fix obs); write conclusions; regenerate SON obs with
  coverage masking for any future multi-year run.
- Forward validation run at the posterior mean vs model defaults (1 GPU run).

## Phase 1 — the identifiability map (the core deliverable, ~2-4 weeks of sweeps)

Build the response matrix **R[parameter family × observable]**: sign, magnitude
vs noise floor, and cross-parameter degeneracy, measured on a fixed platform
(single-year Oct, 7-day spin-up, Atmos 0.42.2, corrected preprocessing).

**Rows (parameter families), in priority order:**
1. **Warm rain** (q_liq_threshold, rain_τ) — ✅ complete (retrospective sweeps).
2. **EDMF entrainment/detrainment** coefficients — NEXT. Hypothesis: moves cl and
   swcre (the "too few" axis) with a different direction from warm rain.
3. **Cloud optics / droplet number / effective radius** (radiation microphysics) —
   hypothesis: moves swcre at ~fixed lwp (the "too bright" axis).
4. **Cloud fraction scheme** (e.g. quadrature/eps params) — targets cl directly.
5. Later: **ice microphysics** (targets lwcre, clivi), **mixing length/diffusion**
   (targets ta/hur profiles — pressure-level obs already supported).

**Columns (observables):** lwp, pr, cl, swcre, lwcre (+ ta/hur@850/500/200 hPa when
mixing-length row is added). Noise floors already measured.

**Method per row (the sweep-first gate, now standard):**
- 3-5 point 1-D sweep per parameter (or small 2-D for tightly coupled pairs) via
  `run_parameter_sweep.jl` — ~5 GPU members × ~2h each.
- Analysis gates (same as the ones pr passed and swcre/cl documented):
  swing/noise per obs, Pearson direction per param, reachability, misfit/noise at
  best. Where diagnostics are already saved from prior runs, prefer **free
  retrospective analysis** over new sweeps (this worked for pr, swcre, lwcre, cl).
- Record each row in a versioned artifact:
  `experiments/calibration/amip/identifiability_map.md` (+ a JLD2 with the numbers),
  stamped with the model version (map rows are Atmos-version-specific; response
  *directions* usually transfer across upgrades, absolute biases do not — re-verify
  rows retrospectively after any model bump).

**Admission criteria distilled from experience:**
- A parameter enters the joint set iff some observable responds with swing/noise
  ≳ 0.3 (lwp's level) and it is not >0.9-degenerate with an existing parameter on
  ALL responsive observables.
- An observable enters iff it is reachable (obs mean inside swept model range) OR
  its structural bias is honestly floored, AND at least one parameter moves it above
  its noise.

## Phase 2 — mid-scale joint calibration (first assembly)

- 4-6 parameters (warm-rain pair + entrainment/detrainment + optics knob) ×
  {lwp, pr, cl, swcre}. TransformUnscented 2p+1 → 9-13 members: still cheap.
- **State falsifiable predictions up front** from the map (as "pr will pull rain_τ
  down" was predicted and confirmed): e.g. "entrainment will be constrained by cl,
  swcre by the optics knob, q_liq unchanged from Phase 0 posterior."
- Run packaged (PACKAGING_PLAN.md) with steering indicators
  (STEERING_INDICATORS_PLAN.md) watching every iteration.
- Success = every parameter's landing point attributable to a mapped response;
  every persistent >2σ residual attributable to known structure.

## Phase 3 — statistical scaling

- Multi-season / multi-year minibatch targets (SON machinery, corrected obs) so the
  larger parameter set is not fit to one month's weather.
- Expand vertical/profile observations (ta/hur pressure levels) with the
  mixing-length family.
- If p grows beyond ~15-20: evaluate EKI variants (localization/sparsity) vs
  TransformUnscented member cost.

## Phase 4 — the large packaged calibration (the boss's objective)

- ~10-20 parameters × ~6-10 observables × multi-year minibatched targets.
- Fully packaged (frozen snapshot + provenance), indicators-instrumented, designed
  directly from the identifiability map.
- The map is what makes it defensible: for every parameter, "this obs constrained
  it, in this direction, at this strength"; for every stubborn residual,
  "structural, here is the missing physics."

## Cross-cutting workstreams (separate plans)

- **Packaging/freeze framework** — PACKAGING_PLAN.md (prevents environment drift;
  becomes mandatory at Phase 2+).
- **Steering indicators** — STEERING_INDICATORS_PLAN.md (automates the per-obs
  fit-to-noise, collapse, degeneracy, reachability verdicts; validated by
  re-discovering this study's manual findings).
- **Model-version discipline**: calibrate only on a frozen model version; after any
  upgrade, re-run the retrospective map analysis on existing outputs before trusting
  old rows (the 0.41→0.42 lesson).
