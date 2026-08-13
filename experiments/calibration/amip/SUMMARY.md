# AMIP calibration campaign — summary, timeline, and scoreboard

A narrative companion to `README.md` (how to launch) and
`CALIBRATION_GUIDE.md` (methodology), written for someone new to this
work. It reconstructs what was tried, in what order, what succeeded, and
what failed.

Sources: the config headers in `config/` (each records that run's design,
falsifiable predictions, and a dated verdict), `CALIBRATION_GUIDE.md`, and
the analysis code, all as of commit `993eec54` ("Add AMIP calibration
machinery", 2026-08-07). Runs whose verdicts are not written in their own
header are reconstructed from cross-references in later headers; the full
run-by-run record (including the sweep driver and `identifiability_map.md`)
lives on the `ne/calibrate` campaign branch.

---

## 1. The problem, briefly

The coupled AMIP model (ClimaAtmos + land, with prescribed sea surface
temperature and sea ice) has cloud and precipitation biases. The campaign
tunes a small set of physical parameters — warm-rain microphysics,
cloud-fraction closure shape, phase-change timescales, terminal
velocities, cloud optics — so that one-to-three-month simulations match
satellite observations: liquid water path (`lwp`, MAC), cloud cover
(`clt`/`cl`, CALIPSO/CloudSat), precipitation (`pr`, GPCP), and cloud
radiative effects (`swcre`/`lwcre`, CERES).

The machinery is Ensemble Kalman Inversion (EKI, specifically
`TransformUnscented` from EnsembleKalmanProcesses.jl). Each **iteration**
runs an **ensemble** of `2p + 1` model realizations (**members**) for `p`
free parameters, compares each member's output to the observation, and
updates the parameter distribution. The **prior** is what you believe
before; the **posterior** is what the data leaves you with. One member is
a ~38-day AMIP run (~1.2 GPU-hours on Derecho), so a typical run of 13–17
members × 6 iterations costs on the order of 100 GPU-hours — which is why
so much of the campaign is about *not* launching runs that cannot work.

Two protocol families appear throughout:

- **Subseasonal**: members start from ERA5 reanalysis initial conditions,
  spin up 7 days, and are graded on the next month(s). The large-scale
  weather is pinned, so parameter signals are cleaner.
- **Free-running AMIP**: members start from a generic initial state and
  run months; the model drifts into its own climate. This is what
  production runs experience, and posteriors do *not* transfer unchanged
  between the two (measured directly; see the ladder in
  `CALIBRATION_GUIDE.md` §8).

Observations are never compared on the raw grid: they are **reduced** to
zonal means (average around each latitude circle) or coarse 2-D blocks,
and the assumed noise includes a **model-error floor**
(`model_error_scale`, a fraction of the field mean) representing the bias
no parameter value can remove.

---

## 2. Read this before trusting any error number: loss ≠ RMSE

Two different "errors" appear in the outputs, and **they routinely
disagree — a falling loss is neither necessary nor sufficient for a
falling RMSE.**

**The EKI loss** (`EKP.get_error`, plotted by `calibration_report.jl`,
"obs RMS x.xx σ" in the steering block) is the misfit the optimizer
actually minimizes: the squared difference between the ensemble-mean
simulation and the observation, computed on the *reduced* (zonal or
block-averaged), *coverage-masked* observation vector, *whitened* by the
assumed noise covariance (interannual variability plus the model-error
floor, sometimes spatially correlated), summed *jointly over all
observables*, against *this iteration's minibatch sample*.

**The ClimaAnalysis-style RMSE** (`post_analyze_iteration.jl` bias
panels, `calibration_report.jl` tables, or any leaderboard comparison) is
a plain cos(latitude)-area-weighted root-mean-square of `sim − obs` in
physical units, per variable.

Documented ways they diverged in this campaign:

1. **The joint loss can fall while one variable's RMSE rises.** The EKI
   trades observables against each other. In the first
   `lwp_clt_swcre_5d` run the update improved lwp and swcre while `clt`'s
   map RMSE *degraded* 0.201 → 0.218 and its mean bias flipped sign
   (+0.078 → −0.117). In `sep_multi`, swcre improved 35 % *paid for* by
   clt degrading 29 %. Loss down, clt worse.
2. **The noise model discounts whole directions of error.** Whitening
   divides each error direction by how much noise you claimed it has. An
   800 km correlated floor treats smooth, large-scale error as "expected",
   so a *near-global* clt bias cost almost nothing in the loss: clt's
   residual was 0.62 σ on the covariance diagonal but only 0.35 σ
   whitened (`clt_only_vterm` verdict). Plain RMSE sees that bias at full
   size; the loss barely does.
3. **The loss target can change every iteration.** With multi-year
   minibatching (`sep_multi`, `phase4`/`phase5`), each iteration grades
   against a *different year*, so consecutive loss values are not the
   same question and small rises/falls between them mean little.
4. **The loss is non-monotone by design.** With the fixed
   `DefaultScheduler(0.1)` step, relax, vterm, and sep_multi all show the
   same shape: loss minimum around iteration 3, rising afterwards. A
   rising tail does not mean the physical fit got worse.
5. **Different spatial samples.** The loss lives on the masked, reduced
   vector (e.g. MAC lwp is ocean-only; masking flips the sign of the lwp
   bias: 11.5 % low masked vs 9.1 % high unmasked). A global unmasked
   RMSE measures a genuinely different quantity.

Practical rule: judge a run by the **per-observable residuals in σ**
(steering block), the **bias maps / g_vs_obs plots** (does the model
track the observed pattern?), and **parameter mean displacement in
prior-σ units** — not by the scalar loss trace, and not by spread
contraction (see the vterm lesson in §4).

---

## 3. The observation vector and the noise covariance: every change, and why

The single biggest evolution in this machinery is *statistical*, not
physical. The original pipeline compared normalized fields on the full
grid with `ScalarCovariance(1.0)` — "every grid point is an independent
measurement with unit noise". That assumption is wrong in both directions
at once, and each change below exists to answer one of two questions
honestly: **how much of the model–obs mismatch is actually noise?** (the
covariance) and **how many independent constraints does the data really
contain?** (the reduction and masks). All of this lives in
`generate_observations.jl`, `preprocessing.jl`, and `correlated_noise.jl`,
mirrored on the simulation side by `observation_map.jl`.

### The covariance (what counts as noise)

1. **`SVDplusDCovariance` replaced `ScalarCovariance(1.0)`.** The noise is
   now *estimated from data*: a low-rank part from the interannual spread
   of the observations across `COVARIANCE_DATE_RANGES` (e.g. October
   2006–2010; rank ≤ number of dates − 1) plus a diagonal part. With the
   scalar unit covariance, EKP had no idea how much mismatch was signal —
   the parameter signal drowned and the error trajectory stayed flat.
2. **Covariance dates were separated from sample dates.**
   `sample_date_ranges` defines *what is calibrated against*;
   `COVARIANCE_DATE_RANGES` defines the realizations the noise is
   estimated from. You can broaden the noise sample (more years) without
   changing the target. Constraint: every sample date must appear among
   the covariance dates (an SVDplusD requirement), which is why the
   earliest usable October is 2006 (CALIPSO coverage starts 2006-06).
3. **A structural model-error floor was added to the diagonal**:
   `(model_error_scale × field mean)²`. This encodes the bias a *perfect*
   parameter set would still leave, so the calibration doesn't distort
   reachable parameters chasing unclosable error. Validated zonal values:
   0.2 for lwp/pr, 0.5 for cl/swcre, 0.25 for lwcre; ~1.5–2× those on 2-D
   grids. Getting this wrong in either direction was a recurring failure:
   too small → collapse (2d10), too generous → the observable goes inert
   and gets traded away (lwcre in phase2 at 0.33 σ, clt at 0.32–0.44 σ in
   phase5 and the first 5d/10d runs).
4. **Per-variable noise groups** (`OBS_NOISE_GROUPS`): each group of
   variables gets its own floor and becomes its own covariance block; the
   per-sample observations are then merged with
   `EKP.combine_observations`. This is how clt got its own floor after
   sharing one with swcre proved 2× too generous. The pipeline's default:
   0.2 for most variables, 0.5 for radiation-like ones.
5. **An optional spatially correlated floor** (`correlated_noise.jl`,
   `decorrelation_length` per noise group): the diagonal floor `D` is
   replaced by `sqrt(D) · K · sqrt(D)` with
   `K_ij = exp(−d_ij / L)` over great-circle distances. The diagonal (per
   -point variance) is unchanged; points in different variables, vertical
   levels, or times stay uncorrelated; 5 % identity is blended in for
   conditioning; the result is stored as a dense `Symmetric` matrix
   (EKP's `SVDplusD` type requires a `Diagonal` floor). The idea: a
   cluster of points inside one L-sized patch counts as ~one constraint.
   The campaign's verdict on it is nuanced — it is a knob with sharp
   edges:
   - at 10° block spacing (1100 km, L = 800 km) it measurably did
     *nothing* (same-seed A/B identical collapse);
   - at 5° spacing (555 km) neighbors correlate at 0.5 and it is kept on;
   - under *zonal* reduction it **double-counts** (averaging already
     removed the small-scale noise) and produced the so_jan null result;
   - on clt it *absorbed a real near-global bias* as "expected noise"
     (0.62 σ diagonal vs 0.35 σ whitened), making clt inert (clt_only).
6. **Per-variable regularization**: `QuantileRegularization(0.05)`
   replaces any global conditioning floor, because variables now sit in
   raw physical units of wildly different magnitude (ta ~ hundreds of K,
   lwp ~ 0.1 kg m⁻²). Knob: `OBS_REGULARIZATION_QUANTILE`; a quantile `q`
   needs ≥ 1/q values per variable, so small reduced vectors (a band's 8
   zonal values) need a larger quantile.
7. **Latitude weighting kept** (`use_latitude_weights`,
   `min_cosd_lat = 0.1`): constraints are weighted by cos(latitude) so
   equal-latitude bands count by their true area, with a floor so the
   poles aren't infinitely discounted.

### The observations (what counts as a constraint)

8. **Normalization was removed on purpose.** Per-variable normalization
   (and its `normalization_stats.jld2`) existed to make the unit
   covariance tolerable. The SVDplusD covariance carries each variable's
   physical scale itself, so normalization is unnecessary and unsupported
   with it; `generate_observations.jl` now actively deletes a stale stats
   file, and the sim side skips normalization when the file is absent.
9. **Spatial reduction is mandatory** (`reduce_spatial`): zonal mean by
   default (NaN-aware `average_lon`); a config-defined `COARSEN_FACTOR`
   switches to cos(lat)-weighted, NaN-aware block averaging
   (`coarsen_lonlat`; never interpolation, which isn't NaN-aware). The
   rationale is quantitative: O(10⁴) correlated grid points against a
   handful of interannual noise modes over-informs a few-parameter
   inverse so badly that the ensemble collapses to a point after a single
   update (spread → 1e-13). Newest option, from the release verdict:
   `ZONAL_SHORT_NAMES` reduces *specific* variables zonally while others
   keep the 2-D grid — used to keep an observable's fixable amplitude
   constraint while deleting an unfixable pattern term (the zsw design).
10. **Missing data is harmonized *before* any reduction**
    (`harmonize_nan_mask_over_dates!`): one shared NaN mask — the union
    of missing points over exactly the covariance dates — is forced onto
    every date. A point missing in any year has no interannual sample and
    must be excluded from all years. Doing this before the zonal mean
    (it used to run after) means every year's zonal mean averages the
    identical set of points, so satellite coverage changes between years
    can no longer masquerade as interannual variability in the estimated
    covariance, and all SVD samples have equal length.
11. **Coverage masks are saved and mirrored onto the simulation**
    (`coverage_mask` → `coverage_masks.jld2` → `apply_coverage_mask` in
    `observation_map.jl`): satellite retrievals are not global (MAC lwp
    is ocean-only — NaN over ~54 % of points), and the mask restricts the
    *simulation* to the same points before *its* reduction. Skipping this
    compares an ocean-only observed mean against an all-longitude
    simulated mean — for lwp that flips the sign of the area-weighted
    bias (11.5 % low vs 9.1 % high). Never bypass it; never impute.
12. **New observables through the loader catalog**
    (`CalibrationTools.CompositeDataLoader`): CALIPSO/CloudSat (`cl`,
    `clt`) and GPCP (`pr`) loaders were added, with explicit
    disambiguation (`lwp` from MAC, not MODIS). `cl` is compared on
    nearest `ALTITUDE_LEVELS` (both sides), `ta`/`hur` on
    `PRESSURE_LEVELS`; a configurable `lat_window()` replaces the
    hardcoded −90…90 (used by the Southern Ocean band runs); and the sim
    side converts units to match the observations (`standardize_sim_units`:
    cl/clt percent → fraction, pr → mm/day).
13. **The observation vector moved into the experiment's `output_dir`**
    (previously saved into the repository), so different setups keep
    independent observations, and the coverage masks travel next to it —
    the pipeline copies them along and warns if they are missing.

The simulation side (`observation_map.jl`) applies the *same* sequence —
level selection, unit standardization, latitude window, coverage mask,
then `reduce_spatial` — in the same order, because the flattened
simulation vector must align positionally with the observation vector,
element by element. Preflight checks much of this (mask applicability,
dimension-name mismatches like CALIPSO `height` vs model `z`, sample
dates present in the data) before any GPU time is spent.

---

## 4. Timeline

Dates are from the config headers and guide. Runs before late July 2026
are undated in this commit; their order is reconstructed from
cross-references.

### Act 1 — LWP-only beginnings and the zonal baseline (early campaign)

The first target was `lwp + cl` (`pressure_levels.jl`, still the
default config in `run_calibration.jl`), motivated by a glaring model
bias: the autoconversion threshold `q_liq` defaulted to 0 — clouds rain
out at the first drop of condensate — leaving LWP ~3× too low.

- **A noise investigation killed `cl` as a target**: its model–obs misfit
  (~0.083) was *smaller* than its own year-to-year variability (~0.10),
  so it carried almost no usable signal. LWP's misfit (~0.031) was well
  above its noise (~0.018). → `lwp_only.jl`. **Succeeded**: parameters
  converged reproducibly, though the error trajectory was
  weather-dominated, with early lwp-only answers around q_liq ≈ 5–7e-4.
- **Longer windows to beat weather noise**: Oct+Nov (`lwp_only_on`),
  Oct–Dec (`lwp_only_ond`, first attempt killed because the prescribed
  SST/sea-ice forcing window ended; fixed by downloading extended forcing,
  `lwp_only_ond_ext`), Sep–Nov (`lwp_only_son`), then multi-year
  minibatching over 2006–2010 (`lwp_only_son_multi`).
- **The spin-up lesson** (`lwp_only_son_multi_spinup`): starting members
  directly from an ERA5 initial condition (which carries no prognostic
  cloud liquid) biased q_liq to ~3.3e-4 vs ~7e-4 with a 7-day spin-up —
  the equilibration transient contaminates the monthly mean. 7-day
  spin-up became standard.
- **Adding precipitation** (`lwp_pr`): q_liq and the rain timescale
  `rain_tau` both raise LWP, so LWP alone cannot separate them; `pr` is
  the autoconversion *sink* and broke the degeneracy. **This became the
  validated zonal baseline: q_liq ≈ 2.8e-4, rain_tau ≈ 1450 s, both
  observables at their noise floors** — the "known answer" later runs use
  as a check.

### Act 2 — going 2-D, and the collapse saga

Zonal means throw away geography (a Pacific bias can cancel an Atlantic
one). The obvious fix — calibrate on a coarse 2-D grid — failed
repeatedly and taught the campaign its noise model:

- `lwp_pr_2d10` / `lwp_pr_2d15` (10° and 15° block averages, diagonal
  noise): **collapsed**. Treating every block as an independent
  constraint over-informs the inverse; parameter spread hit ~0 by
  iteration 5 with residuals still stuck at 2.1–2.4 σ. A collapsed run is
  not just overconfident — it can be *wrong* (the 2d10 run froze q_liq at
  the prior mean).
- `lwp_pr_2d10_corr` (add an 800 km spatially correlated noise floor):
  **still collapsed**, tracking the diagonal run almost exactly — at
  1100 km block spacing the correlation kernel inflates the smooth error
  directions only 2–3×, not the ~10× needed.
- `lwp_pr_2d10_corr04` (prepared 2026-07-29): the honest fix is the floor
  *magnitude* — 2-D residuals stick near twice the zonal floor, so use
  ~2× the zonal `model_error_scale`. This scaling is what later 5° runs
  adopted.

Lesson (now `CALIBRATION_GUIDE.md` §3): start zonal; validate a known
answer before trusting 2-D; floors ~1.5–2× when you go 2-D.

### Act 3 — the phase ladder: more observables, zonal (through 2026-07-31)

A parallel track expanded the *observable* set on safe zonal means,
gated by parameter sweeps (`identifiability_map.md`, on `ne/calibrate`):

- **phase2** (4 parameters × lwp, pr, swcre, lwcre): ran to completion,
  but `lwcre` sat at 0.33 σ from iteration 1 — its assumed floor (0.5)
  was ~2× the real error, so it contributed no constraint.
  → **phase2_lwcre025** halved the floor to 0.25 and validated it.
- **phase3** (6 parameters × 5 observables; adds the cloud-fraction pair
  `eps_rel`/`steepness` and `cl`): **finished healthy 2026-07-30**. Its
  posterior became the campaign anchor: eps_rel 0.0335, detrainment
  coefficient dmfvd ≈ 0.56, detr_buoy 1.61, q_liq in the 2.2–2.5e-4
  zonal cluster.
- **phase4_ond / phase4_son** (3-month members, target year rotating each
  iteration): OND designed around a CALIPSO data gap (no Dec 2009); SON
  staged, then **skipped by decision 2026-07-30** — parked as an ablation.
- **phase5_clt** (phase3's 6 parameters × 6 observables incl. `clt`,
  Sep–Nov members rotating over 2006–2010): **completed 2026-07-31.**
  Passed: no collapse, q_liq 2.33e-4 (zonal cluster confirmed
  seasonally robust), rain_tau finally tightened (0.33× prior).
  Partial: the rotating year kicked parameter drift above threshold at
  iterations 3–4 before settling. Failed: clt fell to 0.32 σ — its zonal
  floor (0.25) still ~2× too generous. New signal: residuals *grow with
  forecast lead time* (pr: Sep 0.72 → Nov 1.43 σ), i.e. the model drifts
  away from reanalysis-pinned states.

### Act 4 — clt on fine grids, the trade-away discovery, and the prior bug (2026-07-28 → 08-02)

- `lwp_cl_swcre_{zonal,2d10,2d15}` tested whether `cl` becomes reachable
  once cloud-fraction parameters are free. The 2d10 diagonal run
  over-informed again (0.01–0.03× spread by iteration 4). A same-seed A/B
  test of the correlated floor at 10° showed it **did nothing**. (The
  zonal companion was healthy but was killed once by a watchdog counting
  only *connected* workers — hence today's long `empty_pool_timeout`.)
- The target switched from `cl` (level-resolved; its CALIPSO definition
  is detection-thresholded and doesn't match the model's quadrature
  cloud fraction) to `clt` (total column cover) on a **5° grid**
  (`lwp_clt_swcre_5d`; 10° merged the coastal stratocumulus decks).
- **The trade-away discovery**: in the first 5d/10d runs clt was 2×
  under-floored, sat at 0.33–0.44 σ, steered nothing — and the optimizer
  *spent* it: clt's map RMSE degraded 0.201 → 0.218 and its bias flipped
  sign while lwp/swcre improved and the joint loss fell. (§2, item 1.)
  → reruns `_cltf035` / `_cltf025` gave clt its own noise group.
- **lwp_clt_swcre_5d_relax** (8 parameters, adding two phase-change
  timescales and the new cloud-fraction floor-release margin; completed
  **2026-07-30**): all four predictions passed. The margin is
  identifiable (the promotion evidence for ClimaParams); the timescales
  landed apart (condensation–evaporation ≈ 106 s,
  sublimation–deposition ≈ 278 s); q_liq 1.93e-4.
- **THE PRIOR BUG (found 2026-07-31)**: EKP's `constrained_gaussian`
  silently returns a *unit normal over the bounds* for small-magnitude
  targets, so **every q_liq prior before 2026-07-31 was the same wide
  prior regardless of what was requested** — prior tightenings never took
  effect. Fix: `checked_constrained_gaussian` (`prior_tools.jl`),
  validated 2026-08-01 in a cheap perfect-model SCM experiment that
  reproduced the bug in miniature. Silver lining: reaching low q_liq
  *against* an effectively wide prior centered at 5e-4 is independent
  evidence the data wants it low.
- **relax_amip** (same problem, free-running AMIP members): measured the
  subseasonal→AMIP shift directly — eps_rel roughly doubles (0.05 →
  ~0.10) and even q_liq moves. This is the basis for the "calibration
  ladder" transfer rules (guide §8): process-local parameters transfer
  with ~1.5× inflated priors; environment-compensating ones (eps_rel, the
  margin, rain/snow timescales) need 3× or wider.
- **The icabl ablation pair** (identical Sep+Oct observations, only the
  initialization mode differs; SUB arm stopped converged **2026-08-02**):
  q_liq 1.34e-4 ≈ relax_amip's 1.35e-4, so the apparent initial-condition
  split in q_liq had been *the broken prior*, not physics — **q_liq is
  IC-insensitive**. eps_rel's shift decomposes ~⅓ season/window, ~⅔
  initialization mode. cond_evap is protocol-invariant (106/117/106 s) —
  the most transferable parameter in the set.

### Act 5 — terminal velocities, the Southern Ocean detour, multi-year September (2026-08-04 → 08-05)

- **vterm** (2026-08-04; three fixed terminal velocities opened, five
  relax parameters frozen): rain fall speed identified (5.42 m/s); snow
  jumped +1.45 prior-σ to 1.69 m/s (see below); ice untouched — the data
  carry no ice-velocity information, so its posterior "is the prior
  talking" and gives **no** evidence for changing the ClimaParams
  baseline. **Method lesson**: spread contracted to 0.04–0.07× prior for
  *every* parameter, informative or not — spread ratio tracks the
  scheduler, not information. **Use mean displacement in prior-σ, not
  spread ratio, as the identifiability test.** Meanwhile clt degraded
  monotonically (0.620 → 0.654) even as total loss beat relax.
- **clt_only** (2026-08-04, **cancelled after 2 iterations** to save ~65
  GPU-hours once the answer was clear): clt got worse *with nothing to
  compete against*, so the joint-run degradation was not competition.
  Root cause: clt already sat well inside its assumed noise (0.35 σ
  whitened vs 0.62 σ diagonal) — the 800 km correlated floor was
  absorbing exactly the smooth, near-global part of clt's error. The
  noise model had declared clt's main bias "expected".
- **so_band** (2026-08-04/05, Southern Ocean 45–65 S band only,
  **stopped after 3 of 6 iterations**): the premise failed — restricting
  to the band made leverage *worse* on all three observables (clt
  residual/spread 5.7 vs 3.7 globally), and v_snow moved −0.24 σ, so the
  global snow signal is *not* a Southern Ocean signal. **Key method
  point: the residual/spread ratio is invariant to the noise scale —
  reweighting cannot buy leverage.**
- **sep_multi** (2026-08-05; five *distinct* Septembers 2006–2010 as a
  minibatch, eps_rel/steepness free again): clt leverage restored (ratio
  2.2; go/no-go passed). Three big verdicts:
  - **v_snow corrected**: +0.60 σ → **1.22 m/s** across five Septembers.
    The October +1.45 σ → 1.69 m/s was substantially one month's weather.
    Quote 1.22, not 1.69.
  - **q_liq replicates again**: 1.329e-4, matching icabl_sub's 1.34e-4 —
    the low-q_liq answer now reproduces across season, years, and three
    parameter sets.
  - **Minibatching does not lower the weather floor** (falsified
    prediction): rotating years removes the *bias* of overfitting one
    realization but does nothing to per-iteration *variance*; that needs
    `minibatch_size > 1` or longer members.
  - **The homogenization warning**: eps_rel moved +2.07 σ (largest of the
    campaign) — *not* because the data wants more saturation variability.
    The optimizer flattened the whole cloud field (model clt collapsed
    from tracking 0.3–0.95 to a flat 0.45–0.6) to dim the too-reflective
    storm tracks: **an optical bias fixed with an areal knob**, swcre
    improving 35 % at clt's expense. Do not quote eps_rel 0.0734 without
    this caveat.

### Act 6 — the optical knob, a null result, and the pattern wall (2026-08-06 → 08-07)

- **so_jan** (2026-08-06; January Southern Ocean, zonal, lwp
  up-weighted): **NULL RESULT — the run learned nothing.** Spread
  contracted 0.7 %; the posterior is the prior. Cause: zonal averaging
  had already removed the small-scale noise, so keeping a 400 km
  *correlated* floor on top **double-counted** the noise and left the
  observation inside its own assumed error (0.31 σ whitened). Fix for any
  zonal rerun: diagonal floor, lower scale. Also flagged: the
  spread-contraction signal-to-noise estimate is an *artifact* unless
  spread actually contracts (< ~0.3) — this run's "S/N" numbers are
  meaningless by construction.
- **release** (2026-08-06/07; frees the four floor-release shape
  parameters plus `Nd`, the prescribed cloud droplet number — the true
  *optical* knob — with clt up-weighted): the designed mechanism worked.
  **Nd moved −0.94 σ to 4.1e7 m⁻³ (41 cm⁻³): dimmer clouds via larger
  droplets instead of fewer clouds.** swcre hit its target (0.553 ≤
  sep_multi's 0.556), eps_rel's excursion halved (+1.17 σ), and
  abs_margin — the parameter added specifically to break the
  eps_rel coupling — moved +0.56 σ, exactly as the closure reading
  predicted. But clt *still* degraded (+22 %, vs sep_multi's +29 %; in
  absolute units 11 % better at equal swcre — the trade got cheaper, not
  free). Posterior leverage: clt ratio 5.6, **out of reach of all six
  parameters**. Conclusion: the remaining clt error is spatial *pattern*,
  and no global scalar in the microphysics or cloud-fraction closure can
  move it — it points at the EDMF (turbulence/convection) state itself.
  Caveat for quoting: 41 cm⁻³ is marine-pristine low; optics may be
  absorbing an EDMF cloud-state error.
- **zsw** (staged, no verdict in this commit): one change from release —
  grade swcre on its *zonal mean* only, deleting the unfixable pattern
  term that finances cloud removal, while lwp and clt keep the full 5°
  field. Reuses release's iteration-1 members (same seed/priors), so one
  of six GPU iterations is free.
- **edmf** (staged, no verdict in this commit): the first parameter set
  that *can* move the clt pattern — six EDMF parameters (detrainment,
  entrainment, static stability, eddy viscosity, max updraft area) mined
  from sweeps, each moving column cover 7.8–10.3 percentage points.
  `mixing_length_diss_coeff` deliberately excluded: swept with exactly
  zero response, suspected not wired.

---

## 5. Scoreboard

### Physical results that held up

| Result | Evidence |
|---|---|
| q_liq (autoconversion threshold) must be nonzero, order 1e-4 — the single most replicated result | Estimates walked down as observables were added and the prior bug was fixed: ~5–7e-4 (lwp-only) → 2.8e-4 (lwp+pr zonal) → 2.2–2.4e-4 (zonal multi-observable) → **1.3–1.5e-4** (5° grid, checked priors; reproduced across season, five years, three parameter sets, and IC modes) |
| cond_evap timescale ≈ 100–117 s, protocol-invariant | relax, relax_amip, icabl — the most transferable parameter found |
| Floor-release margin is identifiable | relax: spread 0.03–0.04× prior with the data moving the mean — the ClimaParams promotion evidence |
| Nd (droplet number) ≈ 41 cm⁻³ fixes swcre brightness through optics, not cloud removal | release; quote with the marine-pristine caveat |
| v_rain ≈ 5.4 m/s identifiable; v_snow ≈ 1.22 m/s (+0.60 σ) | vterm, corrected by sep_multi |
| v_ice unconstrained by lwp/clt/swcre | vterm (mean displacement −0.06 σ) |
| dmfvd ≈ 0.56, detr_buoy ≈ 1.6 stable across phases | phase3 and phase5 |
| eps_rel is environment-compensating: ~0.05 subseasonal vs ~0.10 free AMIP, split ⅓ season / ⅔ IC mode | relax vs relax_amip vs icabl |
| Known unresolved tension: zonal answers (q_liq ~2.3e-4, eps_rel ~0.037) vs 5° answers (q_liq ~1.3–1.9e-4, eps_rel ~0.049–0.058) | phase5 vs relax family |

### Approaches that failed, and what each taught

| Approach | Outcome |
|---|---|
| `cl` (level-resolved cloud fraction) as a target | Dropped: misfit below its own interannual noise, plus a definition mismatch vs the model diagnostic |
| Raw or diagonal-floor 2-D grids (2d10, 2d15, twice) | Ensemble collapse; one run froze q_liq at the prior mean — collapsed ≠ merely overconfident |
| 800 km correlated floor at 10° spacing | Did nothing (same-seed A/B identical); at 5° spacing it helps, under zonal reduction it double-counts (so_jan null) |
| Correlated floor globally on clt | Backfired: declared clt's near-global bias "expected noise", making clt inert (clt_only) |
| Narrowing the region to raise signal (so_band) | Leverage got worse; stopped early |
| Reweighting to make an observable reachable | Impossible — residual/spread is noise-scale-invariant. Weighting only arbitrates *priority* between reachable observables (valid in release, provably useless in vterm) |
| Zero spin-up | Biased q_liq low via the cloud equilibration transient |
| Minibatch over years to lower the weather floor | Removes overfit bias only; the per-iteration floor needs `minibatch_size > 1` or longer members |
| DataMisfitController (adaptive EKI stepping) | Collapsed covariance ~3 orders of magnitude in one step; replaced by fixed Δt = 0.1 |
| Every attempt so far to *improve* clt | Flat or degraded in every completed run — amplitude knobs can't fix a pattern error; the staged EDMF run is the first credible lever |

### Method/infrastructure lessons now baked into the code

- Wiring preflight (a silently-ignored parameter name is a hard failure) —
  a sweep parameter with zero response is presumed not wired.
- `checked_constrained_gaussian` after the silent-prior bug.
- Identifiability = **mean displacement in prior-σ**, not spread ratio.
- Coverage masks applied *before* any reduction (skipping them flips the
  lwp bias sign).
- Correct sample indexing in the forward model (`iter`, not `iter + 1`) so
  member output and observation are the same year under minibatching.
- Steering indicators every iteration; go/no-go leverage gate at
  iteration 1; kill runs early when the answer is clear (clt_only,
  so_band saved ~100 GPU-hours combined).
- Freeze the Manifest mid-study: an Atmos upgrade once changed ocean lwp
  by 27 % at identical parameters.
- Derecho: pack 4 workers/GPU node; bind workers to the HSN for
  login-node drivers; long empty-pool timeout (queued ≠ dead).

---

## 6. Where things stand (as of this commit)

Settled enough to use: the low-q_liq warm-rain answer, cond_evap ≈ 100 s,
the margin's identifiability, the Nd optical pathway, v_rain, and the
transfer rules of the calibration ladder.

Open, in rough priority order:

1. **clt pattern** — the wall every run hit. The staged `edmf` run (with
   the `zsw` loss) is the designed answer; `zsw` alone tests whether
   removing swcre's pattern term stops the cloud-removal trade.
2. **The weather floor** — needs `minibatch_size > 1` (multi-date forward
   model, ~5× workers) or longer members; nothing else measured moves it.
3. **Zonal-vs-2D posterior tension** in q_liq and eps_rel — unresolved;
   the so_jan diagnostic caveat means the one run meant to settle
   reduction questions was uninformative.
4. **clt floor tuning** — 0.25 zonal is still ~2× generous (phase5);
   ~0.12–0.15 is the indicated next value.

---

## 7. Epilogue: how the campaign ended (from `ne/calibrate`, after this commit)

Recovered 2026-08-12 from nefrathenrici's clone
(`/glade/u/home/nefrathe/clima/ClimaCoupler.jl`, branch `ne/calibrate` at
`bfb85870`). The staged runs of §6 ran, and the campaign **closed**:

- **zsw (2026-08-07)** — grading swcre zonally removed about *half* the
  clt cost and stopped the swcre overfit, but the optimizer still sold
  clt for swcre's zonal amplitude "at half the price," paid by eps_rel
  (+0.67 σ) rather than Nd (−0.13 σ, barely used). The pattern term was
  only half the story.
- **edmf (2026-08-08)** — flat. All residual moves inside the weather
  floor; no parameter displaced ≥ 0.5 prior-σ; clt leverage ratio 3.7,
  out of reach. Verdict: *September AMIP at 5° cannot see the EDMF
  state*. With warm rain, cloud-fraction closure, optics, and EDMF all
  exhausted, **the AMIP-rung scalar campaign is converged at this noise
  model**; EDMF calibration moves to the SCM rung.
- **panom (2026-08-08)** — last resort: clt's entire weight on its
  zonal-*anomaly* pattern (ANOM reduction). The pattern residual did not
  move (1.41 → 1.40 σ) while the ensemble collapsed; the frozen
  posterior (eps_rel +2.0 σ, the campaign's largest excursion) was NOT
  adopted. Verdict: *the clt pattern error is directionally orthogonal
  to the entire parameter space — model development or the SCM rung.*
- **Attribution (zero-GPU)** — `clt_regime_attribution.jl` composites
  the clt residual against ERA5 regime predictors (lower-tropospheric
  stability, ω₅₀₀): the residual concentrates in the
  stratocumulus-to-cumulus transition, and the B9 alignment test
  (parameter response maps projected onto the residual in whitened
  space) is the number future scoping should gate on.

Where the artifacts now live for this account: campaign-only docs,
analysis scripts, and all 42 configs with final verdicts are in
`campaign_reference/` (this directory). Analysis artifacts of the four
endgame runs — release, zsw, edmf, panom (reports, EKI files, plots,
logs, observation vectors; member model output excluded) — are under
`/glade/derecho/scratch/kphan/nefrathe_calibration_runs/`. The full set
of earlier run outputs remains on nefrathe's scratch
(`/glade/derecho/scratch/nefrathe/amip_calibration_*_out`), reachable
while it survives the scratch purge cycle.

## 8. Glossary

**Observables** (satellite product in parentheses): `lwp` liquid water
path (MAC, ocean-only); `cl` height-resolved cloud fraction and `clt`
total column cloud cover (CALIPSO/CloudSat); `pr` precipitation (GPCP);
`swcre`/`lwcre` shortwave/longwave cloud radiative effect (CERES); `ta`,
`hur` temperature/humidity (ERA5).

**Parameters** (shorthand → meaning): `q_liq` — cloud-liquid
autoconversion threshold, how much condensate a cloud holds before
raining; `rain_tau`/`snow_tau` — autoconversion timescales; `eps_rel`,
`steepness` — cloud-fraction closure shape (subgrid saturation
variability floor and its sharpness); `margin`/`abs_margin`/`sharpness`/
`residual` — cloud-fraction floor-release shape; `cond_evap`,
`subl_dep` — condensation/evaporation and sublimation/deposition
relaxation timescales; `dmfvd`, `detr_buoy`, `entr_coeff` — EDMF
detrainment/entrainment coefficients; `Nd` — prescribed cloud droplet
number concentration (cloud optics); `v_ice`/`v_rain`/`v_snow` — fixed
terminal fall speeds.

**EKI terms**: *member* — one model run at one parameter sample;
*iteration* — one ensemble + one parameter update; *spread* — ensemble
standard deviation, quoted as a fraction of the prior's; *drift* —
movement of the ensemble mean per iteration; *residual (σ)* — misfit
divided by the assumed per-point noise; *whitened* — transformed by the
full noise covariance; *noise floor / `model_error_scale`* — assumed
irreducible model error as a fraction of the field mean; *leverage
ratio* — iteration-1 residual over ensemble spread of the predicted
observation: how far the target is versus how far the parameters can
reach (≳3 means effectively unreachable); *minibatch* — grading each
iteration against a different observation sample; *collapse* — spread
crashing toward zero while residuals are still large, i.e. false
confidence.
