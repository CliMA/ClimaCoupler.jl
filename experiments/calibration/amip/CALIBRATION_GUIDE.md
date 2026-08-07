# Global calibration guide

How to calibrate the coupled AMIP model against satellite observations.
Update it when a run teaches something new. Companion: the campaign
configs in `config/`, whose headers hold each run's design, predictions,
and verdict.

## The workflow at a glance

1. Pick observables and check their definitions against the model diagnostics.
2. Measure each observable's model error floor.
3. Pick parameters by mechanism and confirm they are wired.
4. Sweep before you calibrate. Admit parameters through the gates below.
5. Reduce the spatial dimension to approximately independent constraints.
6. Design priors from the sweep.
7. Run preflight. Submit only on green.
8. Read the steering indicators every iteration.
9. Record the verdict in the config header.

Cost scale: one member is a 38 day AMIP run, ~1.2 GPU hours. A calibration
runs `(2p + 1) * n_iterations` members. The gates exist to avoid spending
that on a run that cannot work.

## 1. Choosing observations

**Ensure the variable definition matches the model diagnostic.** The model
diagnostic and dataset can measure different things under the same name. For
example, CALIPSO `cl` is detection-thresholded radar-lidar cloud fraction, while
the model's `cl` is any-condensate quadrature cloud fraction. Comparing these
two results in a biased calibration.

**Dataset coverage can be missing or biased.** Retrievals are not global (for
example, MAC lwp is ocean-only). Currently coverage masks are saved with the
observation vector and applied to the simulation before any reduction, so both
sides average over the same points. Never bypass it, and never impute missing
data.

**Use the loader catalog.** `ClimaCoupler.CalibrationTools.CompositeDataLoader`
provides a catalog of preprocessed datasets for use in calibration. If your
desired dataset is not found, add a new DataLoader.

**Measure the model error floor first.** The floor is the bias the model keeps
at the optimal parameter set: structural error that no parameter value can
remove. It is added to the diagonal of the SVDPlusD covariance.

To measure the error, run the model at the best parameters you know, then
compute the model's RMSE against the observation (ideally reduced), and
divide it by the observation magnitude: model_error_scale = ‖y − G(θ*)‖ / ȳ.

Validated zonal values: 0.2 for lwp (MAC) and pr (GPCP), 0.5 for cl (CALIPSO)
and swcre (CERES), 0.25 for lwcre (CERES). 2D grids need about 1.5 to 2x
these values to offset the larger point count.

## 2. Choosing parameters

1. **Start from the closure.** Pick the parameters of the closure that controls
the observable's bias. If you are unfamiliar with the details of the closure,
look through the relevant source code and literature.

2. **Sweep before you calibrate.** Notation: `θ` is the parameter vector,
`θ0` the base point, and `G_v(θ)` the reduced field (zonal or
block-averaged) of observable `v`, with `‖·‖` the cos(lat)-weighted RMS
over its points and an overbar the cos(lat)-weighted mean. A sweep of
parameter `p` runs members at `θ_lo` and `θ_hi` (the ends of `p`'s
plausible range, all other components at `θ0`) plus replicate members
`θ_rep` that jitter every component of `θ0` by 0.1 percent (physically
identical, so their output differences are pure internal variability).

Criteria of useful parameters:

- **Does the parameter have an effect on the observations bigger than
  weather noise?** Signal is the swing across the parameter's plausible
  range, `s_v(p) = ‖G_v(θ_hi) − G_v(θ_lo)‖`. Noise is the replicate
  spread, `n_v = mean_i ‖G_v(θ_rep,i) − G_v(θ0)‖`: the replicates are
  physically identical, so their differences are pure internal
  variability. Require `s_v/n_v ≥ 0.3` for at least one observable.
  Below that, the calibration cannot tell the parameter's effect apart
  from weather, and no prior tuning helps — only more data (longer runs,
  more months) changes the ratio.
- **Does the parameter have a distinguishable effect from every other
  parameter in the set?** Each parameter has a fingerprint: its response
  pattern `r(p) = (G(θ_hi) − G(θ_lo)) / (θ_hi − θ_lo)`, stacked over the
  responsive observables. Require `|corr(r(p), r(q))| ≤ 0.9` (Pearson,
  over the stacked points) against every admitted parameter `q`. Above
  0.9 the two push the observations in the same direction, so the update
  learns only their combination and the pair slides along a valley of
  equally good fits. An observable that responds to `p` but not `q`
  breaks the tie — pr broke the q_liq/rain_tau pair this way.
- **Can the model reach the observed value at all?** Require
  `min_θ Ḡ_v(θ) ≤ ȳ_v ≤ max_θ Ḡ_v(θ)`, where the min and max run
  over the swept members and `ȳ_v` is the observed mean. If violated,
  the structural gap `Δ_v = dist(ȳ_v, [min Ḡ_v, max Ḡ_v])` is
  unclosable by `p`; calibrating against it distorts the reachable
  parameters. Require the floor to cover it (`σ_v ≳ Δ_v`, with `σ_v`
  the per-point floor from `model_error_scale · ȳ_v`) or drop `v`.


3. **Design priors from the sweep.** The sweep shows curvature and asymmetry.

## 3. Building the observation vector

**Never calibrate on the raw comparison grid.** Correlated points counted as
independent constraints over-inform the inverse and collapse the ensemble.

**Choose the reduction deliberately.**

- Zonal mean: the default.
- Coarse 2D blocks: keeps geography. Block-average with
  cos(latitude) weights. Never interpolate as it is not NaN-aware.
- Size blocks at or above the decorrelation scale (measured: ~625 km for
  lwp, ~875 km for pr). Both 10 and 15 degree diagonal-floor runs still
  over-informed; start zonal and validate a known answer before trusting 2D.
- A correlated floor (`decorrelation_length` in `OBS_NOISE_GROUPS`, see
  `correlated_noise.jl`) did NOT help at L = 800 km on the 10 degree grid:
  same-seed A/B runs contracted identically from iteration 2 and collapsed
  the same way. At 1100 km block spacing the kernel inflates the smooth
  directions only 2 to 3x, and the gap to zonal is 10 to 16x. The floor
  magnitude is the lever: 2D residuals stick near 2x the zonal floor, so
  use about twice the zonal `model_error_scale` (test in flight).

**Harmonize missing data before reducing.** `harmonize_nan_mask_over_dates!`
forces one shared NaN mask across the covariance dates. A point missing in any
one date has no interannual sample, so it must be masked out of every date.
Build the union-of-missing-points mask from only the dates in
`COVARIANCE_DATE_RANGES` rather than the entire dataset.

**Match dimension names.** Products spell axes differently (CALIPSO `height`
vs model `z`). Preflight fails loudly on a mask that applies to zero points.

## 4. Running the calibration

**Validated EKP settings:** `TransformUnscented` with `impose_prior = true`,
`2p + 1` members, `DefaultScheduler(0.1)` with `n_iterations ~ 1/dt` (8 to
10), fixed `rng_seed`. Do not use DataMisfitController: no step size knob,
and one over-informative observation collapses the covariance in one step.

**Launch procedure:** generate observations, `preflight.jl <config>` (exit 0
required), qsub. The pipeline (`pipeline.jl` / `run_amip.jl`) automates this
order.

**Derecho:**

- The `main` queue charges whole nodes. Drivers go on CPU nodes
  (`select=1:ncpus=16`, `CLIMACOMMS_DEVICE=CPU`); GPU workers pack 4 per
  node (`workers_per_node = 4`, ~26 GB each), cutting charged GPU hours ~3x.
- Small jobs go to `develop`; sweeps can use `preempt` (0.2x charge).
- Give the worker pool a long empty-pool timeout: the watchdog counts only
  connected workers and kills a healthy driver whose workers are queued.

## 5. Reading a running calibration

The driver logs a steering block every iteration:

- `obs RMS x.xx sigma`: residual over the per-point sigma. Above 2 means
  signal remains; near 1 means fit to the floor (success); well below 1
  means the floor is too generous.
- `spread y.yy x prior`: gentle contraction is health. A few percent of the
  prior with residuals above 2 sigma is the over-informed signature;
  COLLAPSED means stop and raise the floor.
- `drift z.zz x spread/iter`: above ~0.15 means not converged. Large drift
  with tiny spread is over-informed again.
- `reach N% in envelope`: persistently one-sided residuals mean structural
  bias.
- `degen`: the runtime version of the sweep's degeneracy gate.

Done means residuals at their floors, stable spread, low drift. Constrained
parameters should reproduce across setups (q_liq came out near 3e-4 in three
configurations).

## 6. Pitfalls that cost real time

Defenses are built in; do not disable them.

1. Raw-grid observations collapse the ensemble. Reduce first.
2. Unmasked coverage flips biases. Masks apply before any reduction.
3. Not-wired parameters burn whole sweeps. Gate on wiring.
4. Definitional gaps are not parameter errors. Compare definitions first.
5. Environment drift invalidates comparisons: a mid-study Atmos upgrade
   changed ocean lwp by 27 percent at identical parameters. Freeze the
   Manifest while runs are in flight; record versions with results.
6. Resuming over partial members writes wrong-dated diagnostics and crashes
   the observation map. Clean partial members first (preflight checks).
7. Killing a driver orphans its workers. Check `qstat` for `julia-*` after
   any kill.
8. A collapsed run is not just overconfident, it can be wrong: the 2d10 run
   froze q_liq at the prior mean.
9. No mutating copies between run directories without a dry run and a
   `parameters.toml` identity check. Never pipe a transfer through `head`.
10. Home quota kills checkpoints. Output goes to scratch via symlinks.
11. Scripts that parse ARGS must `empty!(ARGS)` before including coupler
    code, or the coupler's ArgParse exits on them.

## 7. Reproducibility practices

- One committed config file per experiment: observables, dates, floors,
  reduction, priors, seed.
- Write falsifiable predictions into the config header before launching.
- Record every verdict in the config header with numbers, run
  directories, and dates. Runs get deleted; the record does not.
- Replicate before believing: a constrained parameter should reproduce under
  a different sample, spin-up, or reduction.

## 8. The calibration ladder (multi-scale strategy)

Subseasonal and free-running AMIP are both targets, and posteriors do not
transfer unchanged between them: the forward map differs, so its optimum
differs. The AMIP-mode relax rerun measured this directly. Use the
hierarchy so each level pays only for what it alone can constrain.

The rungs, cheapest first:

1. **Sweeps** (preempt queue): identifiability gates. Does the parameter
   move the observables, distinguishably from the others? (Section 2.)
2. **SCM** (ClimaAtmos `calibration/experiments/gcm_driven_scm`, plus
   BOMEX/ARM/LES column configs; seconds per member): process parameters
   in controlled regimes — autoconversion, relaxation timescales,
   entrainment/detrainment shape.
3. **Subseasonal global** (~2.5 h members): many parameters in the right
   large-scale state. ERA5 initialization pins the environment, so the
   posterior is close to process-true.
4. **Free AMIP** (~7 h members): few parameters, few iterations. The model
   is graded in its own climate, which is what production runs experience.

Information flows down the ladder as priors: the level-k posterior becomes
the level-k+1 prior, inflated per parameter by its measured cross-level
shift. Parameters split into two classes, and the campaign measured both:

| parameter class | evidence | prior at the next level |
|---|---|---|
| process-local (transfers) | dmfvd 0.56 in phase 3 and phase 5; detr_buoy back to 1.6; cond_evap ~100 s in both relax modes | posterior, ~1.5x inflated |
| environment-compensating | eps_rel 0.05 (subseasonal) vs ~0.10 (AMIP); rain/snow tau swings; q_liq slid below the cluster in AMIP | posterior center, 3x or wider |

The steering drift line is the transfer diagnostic: if the next level
pulls a transferred prior harder than ~2 sigma by iteration 2, the
parameter is environment-compensating. Widen its prior and record it as
non-transferable instead of fighting it.

Two prerequisites. Priors must actually take effect: use
`checked_constrained_gaussian` (prior_tools.jl), since a silent
unit-normal prior discards the lower level's information (this happened
to every q_liq prior before 2026-07-31). And do not fully freeze process
parameters at the AMIP level: the AMIP rerun moved even q_liq, so
"process" parameters absorb some environment bias too. Tight beats fixed.

Validated 2026-08-01 with ClimaAtmos's perfect_scm experiment (origin/main,
TRMM column, minutes per member on CPU): the stock run recovers its known
terminal velocities to ~3%, and a variant adding q_liq (truth 1.8e-4)
reproduces the prior bug in miniature — the plain constructor built a
(5.0e-4, 2.1e-4) prior instead of the requested (3e-4, 1e-4) and the
calibration never left 5e-4, while the checked constructor built the
exact prior and moved toward truth. One cheap perfect-model run per new
parameter set is the pipeline's unit test; run it before spending GPU
hours. Caveat from the same demo: the water-path scalars only weakly
constrain q_liq (2.9e-4 after 5 iterations), so an SCM rung is only as
good as its observation design — profile observations or multi-case
forcing are needed before SCM posteriors can serve as priors.

## 9. Tool reference

| Tool | What it does |
|---|---|
| `preflight.jl <config>` | Pre-submission gate: config, observations, run dir, wiring |
| `pipeline.jl` / `run_amip.jl` | Three-input launcher: YAML + priors + observables |
| `generate_observations.jl` | Observation vector and coverage masks from a config |
| `run_calibration.jl` | Calibration driver (CPU node, packed GPU workers) |
| `steering_indicators.jl` | Per-iteration health block in driver.log |
| `post_analyze_iteration.jl` | Per-iteration plots, bias maps, evolution GIFs |
| `calibration_report.jl` | End-of-run report |
