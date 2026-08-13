# Identifiability Map

**The response matrix R[parameter family × observable]** driving the calibration
expansion (ROADMAP.md Phase 1). Each row reports, per observable: parameter-driven
signal vs the observation's own noise, response direction, degeneracy with existing
parameters, and reachability. Rows are **model-version-stamped** — response
*directions* usually survive model upgrades, absolute biases do not (the Atmos
0.41.3→0.42.2 upgrade changed ocean LWP by 27% at fixed parameters).

**Platform** (fixed across all rows): single-year Oct-2010 target, 7-day spin-up,
h_elem=12, subseasonal/WeatherModel mode, zonal-mean observations with coverage
masking (commit 6ffd3646), SVDplusD-style noise = interannual (Oct 2006–2010) + 20%
model-error floor.

**Model version: ClimaAtmos 0.42.2** for every row through zsw (2026-08-07).
The EDMF calibration row runs on **0.42.3** (post-rebase); treat sigma
comparisons across that boundary as loose. Older runs on 0.41.3 are not
comparable.

## Measured noise floors (zonal-mean, Oct)

| observable | source | noise RMS | units |
|---|---|---|---|
| lwp | MAC (ocean-only, 54.6% masked) | 0.026 | kg/m² |
| pr | GPCP (full coverage) | 0.73 | mm/day |
| swcre | CERES | 9.8 | W/m² |
| lwcre | CERES | 5.1 | W/m² |
| cl | CALIPSO (3 altitude levels) | 0.043 | fraction |

## Admission gates (from the warm-rain study)

- **Parameter enters the joint set** iff some observable responds with
  swing/noise ≳ 0.3 AND it is not >0.9-degenerate with an existing parameter on all
  responsive observables.
- **Observable enters** iff reachable (obs mean inside swept model range) or its
  structural bias is honestly floored, AND ≥1 parameter moves it above its noise.

---

## Row 1 — warm-rain autoconversion (q_liq_threshold, rain_τ) — ✅ COMPLETE

Sources: 40-member retrospective sweep (LWP-only run), iteration-1 sigma-point
design, 20-member scout (lwp+pr run, Atmos 0.42.2), lwp+pr calibration itself.

| observable | swing/noise | direction (RMSE vs param) | reachable? | verdict |
|---|---|---|---|---|
| lwp | 0.33 | q_liq r≈−0.88 (pre-fix run); wants q_liq ≈ 3e-4 (corrected) | yes (fit to ~1σ) | **constrains q_liq** |
| pr | 0.76 | rain_τ↓ improves pr 34% while lwp −1.4% | yes (fit to ~1σ) | **constrains rain_τ** (breaks the q_liq/rain_τ degeneracy lwp alone leaves) |
| swcre | 1.98 | q_liq r = +0.95 (wants q_liq as low as possible) | **NO** — obs −40.2 W/m², model (−68.6, −53.3) | unreachable: model too bright at every setting |
| lwcre | 0.19 | dead | no (but mean only ~4% off) | uninformative for warm rain |
| cl | 0.33 | q_liq r = −0.76 (wants q_liq high) | **NO** — obs 0.170, model (0.139, 0.154) | unreachable: cloud fraction 10–20% too low |

**Calibrated result (lwp+pr joint, Atmos 0.42.2):** q_liq → ~2.9–3.1e-4,
rain_τ → ~1150–1460, both observables at ~1σ by iteration 3, no collapse.
(Historical q_liq≈7e-4 is obsolete: coverage-mask sampling bug + model upgrade.)

**Replication (2026-07-28):** three independent setups on Atmos 0.42.2 agree —
single-Oct lwp+pr (masked obs): ~3.0e-4 · SON×5yr minibatch spinup-0: ~3.2e-4 ·
SON×5yr with 7-day spin-up: ~3.4e-4. The spin-up-vs-no-spin-up difference is inside
ensemble spread → **spin-up does not move the posterior** (it improves absolute fit
~15%); the old single-vs-multi-year "divergence" was entirely the masking artifact +
model version. q_liq ≈ 3e-4 is now the replicated warm-rain result.

**Structural finding ("too few, too bright"):** at the warm-rain optimum, clouds are
~18% too few (cl) yet 13–28 W/m² too reflective (swcre), and the two pull q_liq in
opposite directions. No warm-rain setting reaches either → the next parameter
families must move cloud *fraction* and per-cloud *brightness* independently of LWP.

---

## Row 2 — EDMF family (6 knobs) — 🔄 IN FLIGHT

**Hypothesis:** the EDMF closure's active knobs move cloud fraction (cl) and
brightness (swcre) in directions warm rain cannot — the "too few, too bright" axes.

### Which observables each parameter is compared against

**Every swept parameter is scored against ALL five observables** — that is what
makes a map row (surprises are the point: pr's rain_τ sensitivity was found this
way). The table records the *hypothesized primary* targets and mechanism; the sweep
verdict fills the actual responses.

| parameter (base) | mechanism | primary obs expected | secondary |
|---|---|---|---|
| `entr_coeff` (0.1) | entrainment velocity scale — dilutes updrafts; more entrainment → weaker/moister-mixed updrafts | **cl**, **swcre** | lwp (via cloud water), pr |
| `detr_buoy_coeff` (1.0) | buoyancy-driven detrainment — where updrafts deposit condensate/area | **cl** (vertical structure), **swcre** | lwp |
| `detr_massflux_vertdiv_coeff` (0.3) | detrainment from mass-flux convergence — cloud-top outflow | **cl** at 2–5 km | swcre |
| `EDMF_max_area` (0.7) | updraft area cap — direct ceiling on convective cloudy area | **cl**, **swcre** | pr |
| `mixing_length_eddy_viscosity_coefficient` (0.14) | environment SGS diffusion — stratocumulus/PBL mixing strength | **swcre**, **lwp** (Sc decks dominate both) | cl at 2 km |
| `mixing_length_diss_coeff` (0.22) | TKE dissipation — sets equilibrium TKE, hence mixing | **swcre**, **lwp** | cl |
| — `turb_entr_param_vec` | NOT swept: `a·exp(−b·area)` ≈ 0 at config values [1e-4, 1e4] — effectively inactive (also a vector param) | — | — |
| — `entr_inv_tau`, `detr_coeff`, `detr_vertdiv_coeff` | NOT swept: disabled (=0); enabling changes the closure form — a different experiment | — | — |

lwcre is scored too but expected flat for all six (warm-boundary-layer knobs; the
warm-rain row already showed lwcre nearly unbiased).

### Design

Plus design (one-at-a-time around base), **12 sweep members + 1 center = 13 runs**:

| member | perturbed param | value |
|---|---|---|
| 1, 2 | entr_coeff | 0.03, 0.3 |
| 3, 4 | detr_buoy_coeff | 0.3, 3.0 |
| 5, 6 | detr_massflux_vertdiv_coeff | 0.1, 0.9 |
| 7, 8 | EDMF_max_area | 0.4, 0.9 |
| 9, 10 | mixing_length_eddy_viscosity_coefficient | 0.05, 0.42 |
| 11, 12 | mixing_length_diss_coeff | 0.073, 0.66 |
| 13 | none (center/baseline) | all base values |

Warm rain **pinned at the calibrated center** (q_liq = 2.885e-4, rain_τ = 1463,
lwp+pr iteration_005 member_001) — the base TOML's q_liq = 1e-4 is far from
calibrated and responses must be measured in the calibrated regime.

**Execution:** `run_parameter_sweep.jl` (sweep block carries this design),
`CALIBRATION_CONFIG=config/lwp_pr.jl`, output
`/glade/derecho/scratch/nefrathe/amip_sweep_entr_detr`. 13/13 members completed
2026-07-28 03:44; analysis: `analyze_row2.jl` (log: tmp/analyze_row2.log).

### ✅ Row-2 RESULTS (ΔRMSE vs center member, in units of each obs's noise;
negative = improves the fit)

| perturbation | lwp | pr | cl | swcre | lwcre |
|---|---|---|---|---|---|
| entr_coeff 0.03 | −0.10 | −0.28 | −0.01 | −0.00 | +0.23 |
| entr_coeff 0.3 | +0.23 | +0.38 | +0.09 | −0.15 | −0.30 |
| detr_buoy_coeff 0.3 | +0.01 | +0.06 | −0.04 | **+0.48** | **+0.37** |
| detr_buoy_coeff 3.0 | +0.15 | +0.11 | −0.01 | 0.00 | −0.15 |
| detr_massflux_vertdiv 0.1 | −0.12 | −0.11 | −0.02 | +0.11 | +0.05 |
| detr_massflux_vertdiv 0.9 | **+1.56** | **+0.74** | +0.03 | **−0.62** | **−0.35** |
| EDMF_max_area 0.4 | −0.26 | −0.08 | −0.01 | −0.07 | +0.06 |
| EDMF_max_area 0.9 | −0.09 | −0.06 | −0.00 | +0.06 | +0.12 |
| mixing_len c_m 0.05 | −0.05 | +0.08 | +0.01 | −0.03 | +0.18 |
| mixing_len c_m 0.42 | −0.22 | +0.33 | +0.12 | +0.06 | +0.12 |
| mixing_len diss 0.073 / 0.66 | 0 | 0 | 0 | 0 | 0 |

**Parameter verdicts** (gate: swing/noise ≳ 0.3 on some observable):

- **`detr_massflux_vertdiv_coeff` — ADMIT, strongest in family** (lwp 1.56, pr
  0.74, swcre 0.62, lwcre 0.35). The first knob that **trades LWP against
  brightness**: 0.9 degrades lwp/pr while *improving* swcre/lwcre. Partially
  degenerate with q_liq on lwp but separable via the opposite-signed swcre
  response — the pair {lwp, swcre} identifies both.
- **`detr_buoy_coeff` — ADMIT, near-pure radiation knob** (swcre 0.48, lwcre 0.37
  with lwp/pr ≤ 0.15). Asymmetric: the signal is on the LOW side (0.3 worsens
  radiation); 1.0→3.0 saturates. Orthogonal-to-lwp radiation leverage is exactly
  what the joint set lacked.
- **`entr_coeff` — ADMIT (marginal)** (pr 0.38, lwcre 0.30, lwp 0.23).
- **`mixing_length_eddy_viscosity_coefficient` — ADMIT (marginal)** (pr 0.33,
  lwp 0.22 improving at 0.42; the largest cl response in the family, still only
  0.12).
- **`EDMF_max_area` — REJECT at these ranges** (≤0.26, and cl 0.01 despite being
  the hypothesized cl driver).
- **`mixing_length_diss_coeff` — NOT WIRED in Atmos 0.42**: `c_d` is derived
  (`c_m·c_b/Ri_c`); the standalone TOML entry is ignored (members bit-identical to
  center). Follow-up candidates for the same physics: `mixing_length_static_stab_coeff`,
  `mixing_length_Ri_crit`.

**Observable verdicts updated by this row:**

- **cl is essentially immovable by the entire EDMF family** (max 0.12σ). The "too
  few" bias is not addressable by entrainment/detrainment/max-area/mixing at these
  ranges → next candidates are the cloud-fraction scheme row (row 4), and BEFORE
  using cl as a calibration target, check the CALIPSO-vs-model cl definition
  (lidar-detectable cloud vs model cloud fraction) — the bias may be partly
  definitional.
- **swcre**: the family moves it (0.35–0.62σ) in the right direction but the obs
  stays outside the swept range (−40.2 vs model −59.6…−51.1) → usable Phase-2
  target only with an honest structural floor; closing the bias likely needs
  row 3 (cloud optics).
- **lwcre: promoted** — now reachable (23.2 ∈ (19.9, 25.6)) and responsive
  (0.3–0.37σ). Candidate secondary target.
- **lwp**: the EDMF family swings it across the observed value — reachable both
  ways. Confirms warm rain is not the only LWP control; joint calibration must
  include swcre to break the q_liq/dmfvd degeneracy.

**Hypothesis scorecard** (what the design table predicted vs found): entrainment →
cl **wrong** (cl immovable); detrainment → cl **wrong for cl, right for
radiation** (both detr knobs are the family's radiation levers); mixing length →
swcre/lwp **half right** (lwp/pr yes, swcre no); max-area → cl/swcre **wrong**
(inert). The one-at-a-time design made every surprise attributable.

---

## Row 2b — TKE-dissipation controls + combined detrainment — ✅ COMPLETE

6 members (2026-07-28): one-at-a-time `mixing_length_static_stab_coeff`
(0.13/1.2) and `mixing_length_Ri_crit` (0.1/0.5), one combined member
`detr_massflux_vertdiv_coeff` = 0.6, one center.

| perturbation | lwp | pr | cl | swcre | lwcre |
|---|---|---|---|---|---|
| static_stab 0.13 | −0.07 | −0.21 | −0.04 | +0.16 | +0.20 |
| static_stab 1.2 | −0.15 | +0.06 | +0.01 | −0.01 | +0.01 |
| Ri_crit 0.1 | +0.12 | +0.08 | +0.02 | +0.01 | +0.08 |
| Ri_crit 0.5 | −0.17 | −0.06 | +0.03 | +0.12 | +0.22 |
| **dmfvd 0.6** | **+0.16** | +0.01 | −0.02 | **−0.34** | −0.13 |

Verdicts:
- **static_stab_coeff, Ri_crit — REJECT** (max swings 0.21σ / 0.22σ). These do
  change the simulation, so the dissipation channel is wired, but it is a weak
  lever at these ranges.
- **dmfvd response is strongly nonlinear**: 0.6 gives −0.34σ swcre at +0.16σ lwp
  cost, while 0.9 gave −0.62σ swcre at +1.56σ lwp cost. Roughly half the
  radiation gain at a tenth of the lwp damage. **Phase-2 prior for dmfvd:
  center near 0.3–0.4 with an upper bound near 0.7.**
- **Cross-sweep determinism confirmed**: the center member reproduces the Row-2
  center to 4 significant figures on lwp/pr/swcre/lwcre. The cl difference
  (0.1523 vs 0.1385) measures the imputation fix: with fabricated points
  removed, the observed cl mean is 0.1788 and the model deficit is ~22%.

## cl definition check (2026-07-28) — the "too few" bias is REAL

Question: is the ~18% cl deficit physics or a definition mismatch between the
observation and the model diagnostic?

**Observation** (`calipso_cloudsat` artifact = CloudSat+CALIPSO radar-lidar
combined product, monthly 2.5°, 240 m height bins): cloud fraction = fraction of
radar-lidar samples where cloud is *detected* (radar reflectivity or lidar
scattering-ratio thresholds), nadir-only ~1:30 LT sun-synchronous sampling. The
loader uses the "all cases" operations slice (full radar ops through 2006–2010).

**Model** (`cl` diagnostic): quadrature SGS cloud fraction —
`cache.precomputed.ᶜcloud_fraction`, the fraction of the subgrid distribution
holding *any* condensate. No detection threshold, full diurnal sampling.

**Directional verdict**: the definitional gap goes the WRONG way to explain the
bias. Detection-thresholded obs *miss* thin cloud, while the model definition
*counts* any condensate — so model cl is definitionally an overcount relative to
"detectable" cloud. The model sits ~18% BELOW the obs despite that overcounting
definition → **the "too few clouds" bias is real and, if anything, understated.**
(Secondary terms are small at our levels: radar ground clutter affects <1 km,
below our 2 km level; diurnal sampling bias for marine Sc is a few percent.)

**Bug found during the check**: the Calipso loader **imputes missing/NaN obs
values with the global mean** (fabricating data in unobserved bins — polar caps
beyond the ~82° orbit limit and gappy cells) which also blinds the coverage-mask
machinery (cl showed 0% missing in scouts because imputation removed the NaNs
first). Fixed: NaNs are now retained so `harmonize_nan_mask_over_dates!` +
`coverage_mask` handle them like MAC lwp's gaps.

**Consequence**: cl stays a legitimate future target (with an honest floor), but
no swept parameter moves it — the deficit points at the cloud-fraction closure
itself → row 4 (quadrature/SGS-variance parameters) is the mechanism to test.

## Rows 3+ — planned (ROADMAP.md priority order)

3. **Cloud optics / droplet number / effective radius** — targets swcre's "too
   bright" axis at ~fixed lwp.
4. **Cloud fraction scheme** (quadrature/eps params) — targets cl directly.
5. **Ice microphysics** — targets lwcre/clivi. **Mixing length / diffusion** —
   targets ta/hur pressure-level profiles.

---

## Spatial decorrelation of the covariance samples (measured 2026-07-28)

Pair correlations of the interannual Oct 2006-2010 anomalies across the 5 years,
pooled per 250 km distance bin on the 2.5 degree grid (script:
tmp/decorrelation.jl): **e-folding ~625 km for lwp, ~875 km for pr** (r at
1000-1250 km is 0.09 and 0.05). Consequence for 2-D coarse observations:
10 degree blocks (~1100 km) are effectively independent pixels, 15 degree
blocks fully independent, and 5 degree blocks would be 30-50% correlated with
their neighbors. This is the empirical footing for COARSEN_FACTOR = 4 (10 deg)
as the preferred 2-D grid, pending the validation runs.

## Six-run comparison verdicts (2026-07-29)

All six overnight calibrations finished. Final states, iteration 6 (phase2:
iteration 8):

| run | q_liq | rain_τ | residuals | final spread (×prior) | verdict |
|---|---|---|---|---|---|
| phase2 (zonal, 4p×4obs) | 2.32e-4 | 1083 | lwp/pr/swcre at floor, lwcre 0.33σ | 0.22-0.75 | healthy; prediction 6 failed (lwcre floor too generous) |
| lwp_pr_2d10 (diag floor) | 5.41e-4 | 1115 | 2.2-2.3σ, stuck | 0.0 COLLAPSED | over-informed; answer wrong (q_liq froze near prior mean) |
| lwp_pr_2d15 (diag floor) | 3.94e-4 | 1295 | ~2σ, drift 2.5-4.3×spread | 0.0 COLLAPSED | over-informed; milder but also collapsed |
| lwp_cl_swcre_zonal | 2.47e-4 | 1776 | all at floor (0.78-0.99σ) | 0.29-0.87 | healthy; cl is reachable with cloud-fraction params |
| lwp_cl_swcre_2d10 (diag) | 2.19e-4 | 1697 | all at floor | 0.0 | right answer, fake confidence |
| lwp_cl_swcre_2d15 (diag) | 2.52e-4 | 1259 | all at floor | 0.02-0.06 | right answer, overconfident |

- **Hypothesis A confirmed**: phase2's lwcre stayed at 0.33-0.37σ for all 8
  iterations under the 0.5 floor. Follow-up launched: `phase2_lwcre025.jl`
  (floor 0.25, predictions in the header).
  **phase2_lwcre025 verdict (finished 2026-07-29 15:20 MDT, job 6946559): all
  three predictions passed.** (1) lwcre landed at 0.61σ, inside the predicted
  0.6-1.1σ band. (2) Detrainment spreads tightened vs the first run:
  dmfvd 0.55 vs 0.58, detr_buoy 0.71 vs 0.73 ×prior — real but marginal, so
  the extra lwcre information is small. (3) lwp 1.03σ / pr 0.93σ /
  swcre 0.81σ at their floors, q_liq 2.33e-4 (in the 2.2-2.5e-4 cluster),
  rain_τ 1283, no COLLAPSED flag, all pairs separable, drift ≤ 0.1.
  Conclusion: 0.25 is the right lwcre floor; phase3 inherits it.
- **Hypothesis B confirmed**: both diagonal-floor lwp_pr 2-D runs hit the
  COLLAPSED steering flag; the cls 2-D runs reached the zonal answer but with
  ~0 spread (same mechanism, masked because their residuals can reach the
  floor). Root cause: the diagonal floor counts every block as independent.
  Fix implemented: correlated noise floor (`correlated_noise.jl`, L = 800 km),
  follow-ups `lwp_pr_2d10_corr.jl` and `lwp_cl_swcre_2d10_corr.jl` launched.
  **Correlation-only verdict (same day): insufficient for lwp+pr.** The
  correlated run hit 0.0x prior spread by iteration 3 with residuals at
  2.2-2.9 sigma, tracking the diagonal trajectory. Two causes: (1) parameter
  sensitivities are smooth and the kernel only inflates smooth directions
  2-3x at 1100 km block spacing; (2) the stuck 2 sigma residual means the
  honest 2-D structural floor is ~2x the zonal 0.2. The run was killed at
  iteration 3 and replaced in its slot by `lwp_pr_2d10_corr04.jl`
  (correlation + floor 0.4). The cls correlated run stays up because its
  residuals sit at the floor; it tests correlation against overconfidence.
  **corr04 verdict (floor 0.4 + correlation, killed at iteration 5,
  12:45 MDT): the doubled floor softens but does not stop the
  over-informing.** Spread 0.27 -> 0.09 -> 0.03 -> ~0.02x prior (prediction
  "final spread above 0.05" failed); pr reached its floor (1.25 sigma), lwp
  stalled near 1.7 sigma. The frozen answer is q_liq = 6.7e-4,
  rain_tau = 1415 — q_liq lands far from the healthy cluster (2-3e-4), so
  at 10 degrees the lwp+pr pair collapses onto a WRONG answer even at
  floor 0.4. Conclusion: lwp+pr does not calibrate on a 10 degree grid
  with these tools; the observables that work in 2-D are the ones whose
  residuals can reach the floor (the cls family). The clt 5/10 pair runs
  with floors 1.5x higher again and its kill rule armed.
  **cls correlation-only verdict (finished 09:48 MDT): also insufficient.**
  Final spread 0.0-0.01x prior (prediction 1 failed; the diagonal twin
  ended the same), residuals at the floor (prediction 3 held). Posterior
  vs zonal (prediction 2, mostly held): q_liq 2.34e-4 vs 2.47e-4,
  steepness 0.58 vs 0.62, rain_tau 1554 vs 1776 and snow_tau 1359 vs 1551
  (both inside zonal's wide spreads); eps_rel 0.057 vs 0.042 is ~1.3x the
  zonal final spread, the one marginal disagreement. Net: at 10 degrees an
  L = 800 km exponential kernel changes neither failure mode; floor
  magnitude is the dominant lever (being tested by lwp_pr_2d10_corr04).
- **Hypothesis C answered**: in the zonal cls run the cloud-fraction
  parameters moved cl to its floor (0.81σ) and stayed constrained
  (eps_rel 0.58×prior) — cl is a usable observable once
  cloud_fraction_eps_rel and cloud_fraction_steepness_scale are in the
  parameter set. Posterior: eps_rel 0.042, steepness 0.62.
- **Healthy-run consensus**: q_liq 2.2-2.5e-4 across phase2 and all cls runs.
  Note this sits below the earlier zonal lwp+pr answer (2.8e-4); rain_τ is
  only weakly identified everywhere (final spreads 0.75-0.87×prior).

## Row 3 — cloud optics — ✅ COMPLETE (2026-07-30, 06:48 MDT)

7 members (job 6956853 + workers): center = phase3 posterior, Nd
(prescribed_cloud_droplet_number_concentration) at 3e7/3e8, ice re
(ice_cloud_effective_radius) at 15e-6/40e-6, 2 jittered replicates.
Scored vs {lwp, pr, swcre, lwcre, cl, clt} zonal Oct 2010
(config/row3_optics.jl). Caveat: noise floors from n=2 replicates —
s/n values are crude, but the verdicts clear the gate by ~10x.

| obs | signal (RMSE range) | noise (rep std) | s/n | reachable |
|---|---|---|---|---|
| lwp | 0.0047 | 0.0019 | 2.4 | no (obs 0.0864 vs 0.0734-0.0759) |
| pr | 0.12 | 0.053 | 2.2 | no (marginal) |
| swcre | 7.31 W m^-2 | 0.54 | 13.6 | NEARLY (obs -40.2 vs -51.4..-42.5) |
| lwcre | 1.52 W m^-2 | 0.31 | 4.8 | yes |
| cl | 0.0038 | 0.0011 | 3.5 | no (definitional gap) |
| clt | ~0.009 | ~0.0008 | 11.6 | — |

- **Nd — ADMIT.** The Twomey fingerprint as predicted: swcre RMSE 21 at
  3e7 -> 16 at 3e8 (center 17.6) with lwp essentially untouched. The
  orthogonal-to-lwp radiation lever Row 2 asked for.
- **ice re — ADMIT.** lwcre RMSE 7.45 (15um) -> 5.9 (40um), and a strong
  swcre response in the same direction: larger crystals improve both.
  De-loads the detrainment pair from radiation duty. Prior should lean
  larger than the 25um default.
- **swcre floor verdict changed**: Row 2's reachability gap (obs -40.2 vs
  model -59.6..-51.1) justified the 0.5 floor; with optics the attainable
  range reaches -42.5, a 2.3 W m^-2 gap. **Phase6 can lower the swcre
  floor toward ~0.3.**
- Surprise: clt responds to optics (s/n 11.6, small absolute) — the
  McICA cover diagnostic is not optics-blind.
- Center member: the phase3 posterior's global lwp is ~14% below obs
  (0.0743 vs 0.0864), consistent with its 1.04 sigma floor residual.
- First clt noise estimate: replicate RMSE std ~0.0008 (n=2).

## clt-floor rerun verdict (2026-07-30, 05:20 MDT)

Jobs 6955107 (5d, clt floor 0.35) / 6955108 (10d, clt floor 0.25), both
exited 0 after 6 iterations. Same configs as the first pair except clt in
its own noise group.

Final means (5d / 10d): q_liq 1.77e-4 / 1.62e-4, rain_τ 1385 / 1661,
snow_τ 1544 / 1400, eps_rel 0.0577 / 0.0586, steepness 0.597 / 0.624.
Final residuals (5d / 10d): lwp 0.67σ / 0.88σ, clt 0.61σ / 0.75σ,
swcre 0.67σ / 0.81σ. Final spreads: 5d 0.02-0.03×prior (frozen, means
only), 10d 0.08-0.16.

- **10d: all predictions passed.** clt landed at 0.75σ (0.7-1.3 band),
  steered the fit, and the cloud-fraction spreads tightened vs the first
  run (eps_rel 0.13 vs 0.16, steepness 0.15 vs 0.18). Not traded away:
  clt RMSE 0.177 -> 0.181 (+2%, vs +8.5% degradation in the
  generous-floor run) while lwp RMSE halved (0.067 -> 0.035) and swcre's
  global bias went -32.4 -> -1.7 W m^-2. **clt floor 0.25 validated at
  10 degrees; adopted for phase5's zonal clt.**
- 5d: clt 0.61σ, just under the band — 0.35 is still ~15% generous at
  5 degrees (prediction 1 marginal fail there).
- **q_liq arbitration: the low 2D answer is real.** With clt pulling,
  q_liq moved to 1.62-1.77e-4, below the first pair (1.75-1.87e-4) and
  well below the zonal cluster (2.2-2.55e-4). The 2D configuration
  systematically prefers low q_liq; phase5 prediction 2 carries the
  cross-check.
- **eps_rel: both grids converge at 0.058**, vs phase3's 0.0335 — the
  eps_rel tension is now two-sided and phase5 prediction 3 arbitrates.

## Phase 3 verdict (2026-07-30, 02:10 MDT)

Job 6950694, exited 0 after 8 iterations. The merged 6-parameter x
5-observable zonal calibration is the campaign's reference posterior.

Final means: q_liq 2.55e-4, rain_τ 1261, dmfvd 0.536, detr_buoy 1.61,
eps_rel 0.0335, steepness 0.585. Final spreads (×prior): q_liq 0.20,
rain_τ 0.74, dmfvd 0.57, detr_buoy 0.63, eps_rel 0.45, steepness 0.66.
Final residuals: lwp 1.04σ, pr 0.94σ, swcre 0.71σ, lwcre 0.64σ, cl 0.79σ.

- Predictions 1, 3, 4 passed: floors by iteration 2, all pairs separable
  throughout, detrainment spreads tightened vs phase2 (0.57/0.63 vs
  0.58/0.73).
- Prediction 2: 4 of 5. q_liq 2.55e-4 extends the healthy cluster's top
  (2.2-2.55e-4 now). **eps_rel missed low: 0.0335 vs the predicted
  0.04-0.06** (cls-family runs gave 0.042-0.067). The joint fit trades
  eps_rel down when the detrainment pair is present; either the
  cls-family value absorbed detrainment effects, or a mild cross-set
  degeneracy below the 0.9 gate. Watch whether phase5 reproduces 0.033.
- rain_τ stays weakly identified (0.74×prior) even with 5 observables —
  the case for the seasonal minibatch (phase5).

## clt 5/10 degree pair verdict (2026-07-29, 21:16 MDT)

Jobs 6949238 (5d) and 6949239 (10d), both exited 0 after 6 iterations.
First calibrations with clt (CALIPSO column cover vs model McICA clt) and
the first direct test of the floor-scaling rule (5d floors = 1.5x the 10d
values against 4x the points).

Final means (5d / 10d): q_liq 1.75e-4 / 1.87e-4, rain_τ 1769 / 1562,
snow_τ 1374 / 1710, eps_rel 0.067 / 0.060, steepness 0.60 / 0.47.
Final residuals (5d / 10d): lwp 0.66σ / 0.88σ, clt 0.33σ / 0.44σ,
swcre 0.64σ / 0.77σ. Final spreads: 5d 0.03-0.05×prior, 10d 0.10-0.20.

- **Floor scaling validated, slightly overdamped.** Neither run collapsed
  or froze above 2σ (the lwp_pr failure mode); the kill rule never armed.
  The 5d run's residuals fell BELOW floor by iteration 2, so 1.5x was more
  than 4x points needed; its spreads still ended at 0.03-0.05×prior
  (prediction 1 failed at 5d, passed at 10d). 10d is the trustworthy arm.
- **clt floors ~2x too generous in both** (residuals 0.33σ/0.44σ from
  iteration 1): clt steered little; lwp and swcre did the constraining.
  The 5d bias maps show the cost: clt RMSE worsened 0.201 -> 0.218 and its
  global bias flipped from +0.078 (over-cover) to -0.117 (under-cover)
  while lwp RMSE halved (0.072 -> 0.038) and swcre's -33 W m^-2 bias
  vanished (-0.9). An under-floored observable is not neutral; the update
  trades it away. A clt-informative rerun needs floors near 0.35 (5d) /
  0.25 (10d).
- **q_liq lands LOW: 1.75-1.87e-4, below the 2.2-2.5e-4 zonal cluster.**
  The two grids agree with each other within 7%, so the shift is
  systematic, not noise: 2D lwp weighting (or clt replacing cl) pulls
  q_liq down ~20%. Open question for phase3 to arbitrate.
- eps_rel 0.060-0.067 vs 0.042-0.057 in the cls family: same direction,
  modestly higher with clt; steepness 0.47-0.60 brackets the zonal 0.62.
  snow_τ and rain_τ remain weakly identified.

## September campaign (2026-07-31 → 08-07): settled parameters

Runs: phase5 (zonal reference), relax + relax_amip + icabl (floor-release
family, IC ablation), vterm (fall speeds), clt_only, so_band, sep_multi,
so_jan, release, zsw. Full verdicts in each config header; distilled here.

**Settled microphysics** (frozen in `toml/calibration_release_fixed.toml`
and carried into `toml/calibration_edmf_fixed.toml`):

| parameter | value | evidence |
|---|---|---|
| q_liq threshold | 1.35e-4 | icabl triangulation: 1.34 (ERA5) vs 1.35 (default IC) with matched priors — **IC-insensitive**; the earlier 1.9-vs-1.35 "IC split" was the broken wide prior, and the 2.2-2.5e-4 zonal cluster belongs to the pre-clt observable set |
| rain_τ | 1800 | never moves (+0.04σ in sep_multi) |
| snow_τ | 1360 | stable since relax |
| cond-evap τ | 100 | relax family, ~106-117 s replicated |
| subl-dep τ | 500 | relax family (333-608 range, weakly identified) |
| steepness | 0.55 | fixed ON PURPOSE: degenerate with eps_rel, leaks into rain via σ_S_eff |
| v_ice | 0.1 | ignored by the data in four runs |
| v_rain | 5.4 | replicated 5.42-5.44, spread 0.055x prior |
| v_snow | 1.22 | sep_multi CORRECTED vterm's 1.69: the +1.45σ pull was ONE September's weather (+0.60σ over five) |

**Method findings** (cost real GPU time to learn):

- **Weather floor**: a single-month target moves parameters up to ~1.4σ on
  weather alone (v_snow). Minibatching over 5 years does NOT lower this
  floor (sep_multi P-fail); it only averages the verdict. A weather-floor
  estimate is valid only when the ensemble contracts (so_jan guard:
  parameter contraction 14x with G contraction 1.76x in release = valid).
- **Spread ratio is NOT an identifiability test** (vterm): every parameter
  contracts to ~0.05x prior because contraction tracks the scheduler's
  pseudo-time. Use mean displacement in prior sigma instead.
- **Single-observable runs cannot steer** (clt_only: 0.62σ flat, worse
  with nothing to compete against). **Regional narrowing hurts leverage**
  (so_band: every ratio worse than global). **An observation already fit
  under its own noise carries zero information** (so_jan null: RMS 0.31σ
  at iteration 1, posterior = prior).

## Row 4 — cloud-fraction floor-release + optics in calibration — ✅ COMPLETE

The release run (lwp+clt+swcre full-field, 6 params: eps_rel, 3
floor-release shape params, floor_residual, Nd) and zsw (identical except
swcre graded zonally, diagonal floor, scale 0.3) close rows 3-4:

- **floor_release_margin — identifiable** (relax: spread 0.03-0.04x with
  a data-driven mean move; the ClimaParams promotion evidence).
- **Nd — the optical pathway is real** (release: −0.94σ, second-largest
  move, dimmer clouds through more droplets at fixed lwp). But zsw shows
  it is only recruited when swcre PATTERN is in the loss: under
  amplitude-only zonal swcre, Nd stays flat (−0.13σ).
- **eps_rel consolidated**: phase3 0.0335 (zonal) → 0.058-0.067 (2D/clt
  runs) → icabl attribution: 1/3 the gap is the Sep-vs-Oct window, 2/3 is
  initialization → release 0.0583 (+1.17σ, chasing swcre pattern) → zsw
  0.048 (+0.67σ, pattern term deleted). **Current best: 0.048.**
- **clt is OUT OF REACH of every closure+optics scalar** at the posterior
  (release posterior leverage: ratio 5.6, response spread 0.168σ, with
  all six free and clt up-weighted). Its remaining error is spatial
  pattern; global scalars fix amplitude only.
- **The clt-for-swcre trade halves but survives amplitude-only grading**
  (zsw, same-Sep-2006 endpoints): clt +0.08σ vs release's +0.17σ, swcre
  stops at its floor instead of overfitting, the surviving half rides
  eps_rel. Any swcre constraint next to clt still buys reflection with
  cloud amount; the EDMF row tests whether pattern-capable parameters
  break this.

## Row 2c — EDMF cover responses, mined (2026-08-07) — REVISES ROW 2

Row 2 scored the EDMF family against `cl` (CALIPSO radar-lidar,
detection-thresholded, 3 altitude levels) and found it immovable
(≤0.12σ). Mining the same sweep outputs against **column cover**
(max-overlap of the monthly 3-D cloud fraction, validated as a clt proxy
at response correlation 0.81 on the release ensemble; random overlap
useless at −0.01) overturns that for clt:

| parameter | cover response (%RMS at sweep range) |
|---|---|
| detr_massflux_vertdiv | 10.3 |
| static_stab | 8.6 |
| Ri_crit | 8.5 |
| entr_coeff | 8.4 |
| eddy_viscosity | 8.3 |
| detr_buoy | 7.8 |
| EDMF_max_area | 7.8 |
| (anchor: eps_rel) | 12.4 (direct clt 17.9) |

Every EDMF knob moves column cover at 60-80% of eps_rel's strength —
**clt is reachable through EDMF**. The row-2 "cl immovable" verdict was
about the thresholded 3-level observable, not cloud amount. diss_coeff
remains dead (zero response, confirms not-wired). Ri_crit duplicates
static_stab and stays fixed.

## Row 2d — EDMF in AMIP calibration (2026-08-08) — ❌ AMIP-INVISIBLE

The EDMF calibration (`config/lwp_clt_swcre_edmf.jl`, 6 params vs the
zsw loss, ClimaAtmos 0.42.3) ran 6 clean iterations and moved NOTHING:
residuals flat inside the weather floor (clt 0.84 -> 0.86, same-Sep
endpoints), displacements all under 0.5 prior sigma (detr_buoy +0.48,
static_stab +0.39, rest <= 0.27).

Iteration-1 leverage (per-observable response spread / residual):

| obs | residual | response spread | ratio | verdict |
|---|---|---|---|---|
| lwp | 0.65 | 0.35 | 1.8 | reachable, already at floor |
| clt | 0.84 | 0.23 | 3.7 | OUT OF REACH (release scalars: 2.6) |
| swcre | 0.76 | 0.22 | 3.4 | hard to reach |

Row 2c's mining was right in direction (dmfvd is the family's strongest
clt lever, 0.41 sigma response) and wrong in sufficiency: sweep-range
responses do not translate into prior-range reach under this noise
model. September AMIP at 5 degrees cannot see the EDMF state.

**Campaign-level conclusion:** warm rain, cloud closure, optics, and
EDMF are now all exhausted against clt, which sits AT its assumed
floor with ratio > 3 for every family. The AMIP-rung scalar campaign
is CONVERGED at this noise model. EDMF parameters move to the SCM rung
(profile observations, per-regime forcing); the AMIP rung re-enters
with a sharper observable or a floor lowered by a better simulator.

## Row 5 — clt pattern anomaly (2026-08-08) — ❌ CLOSED, CAMPAIGN OVER

Metric scoping over existing ensembles found the clt zonal anomaly at
leverage ratio ~1 (vs 1.3-3.7 full field) and the panom run tested it:
clt graded ONLY on its pattern component, closure/optics free, preempt
workers, ClimaAtmos 0.42.3. Result: the pattern residual did not move
(1.41 -> 1.40 sigma, same-Sep endpoints) while eps_rel ran +2.0 prior
sigma and the ensemble collapsed to 0.0x spread (the anomaly's
near-zero per-point floors over-inform; LESSONS B10). The scoping ratio
was quadrature optimism: response MAGNITUDE reaches the residual,
response DIRECTION does not (LESSONS B9).

**Final campaign statement on clt:** four parameter families against
the full field, EDMF against the same, and closure/optics against the
isolated pattern — the clt pattern error is directionally orthogonal to
the entire available parameter space. No representation, weighting, or
floor fixes that. The fix is model physics (or profile-level SCM
constraints), after which calibration re-enters.

## Residual attribution and alignment (2026-08-10, zero GPU)

Regime attribution of the posterior-of-record clt residual (edmf
iteration-1 center vs the 5-September CALIPSO mean; ERA5 LTS and
omega500 predictors; MAC-derived ocean mask; script
clt_regime_attribution.jl in the session scratchpad):

- Global deficit: cos-mean r = -0.105 clt fraction against an observed
  0.67 (model ~16% short of observed cover).
- **The pattern error is the stratocumulus-to-cumulus transition.** The
  top LTS decile (> 17.8 K, the Sc decks) carries 22% of the squared
  residual (28% of the pattern component) at mean r = -0.20, twice the
  global deficit. Across the tropical omega500 gradient the sign FLIPS:
  relatively too much cloud in weak subsidence (+0.06..+0.10), too
  little in the strongest subsidence (-0.04..-0.08). The model breaks
  the decks up too early and moves the cloud into the trades.
- Feasibility of regime-dependent parameters: LTS + omega500 explain
  26.5% of the pattern variance (2 coefficients). Real but partial; an
  eps_rel(LTS) or entrainment(omega500) form has about a quarter of the
  pattern in reach. The remaining 3/4 is not regime-linear.
- **B9 alignments** (cos angle of signed response with the residual,
  noise-whitened): against the FULL residual, eps_rel 0.755 and
  abs_margin 0.618 — their directions ARE the amplitude deficit, which
  explains every eps_rel excursion in the campaign. Against the PATTERN
  component, every parameter in both families is <= 0.27 (best:
  abs_margin 0.273, eddy_visc 0.203, eps_rel 0.195). The best single
  direction can remove at most cos^2 ~ 7% of the pattern misfit, which
  is why panom was flat. Gate future runs on alignment >= ~0.5.

Consequences: SCM cases should target the Sc and Sc-to-Cu transition
regimes (strong-inversion subsidence columns); the model-development
target is deck persistence under strong inversions (cloud-top
entrainment and mixing at the coarse vertical grid); COSP-consistent
diagnostics would sharpen the same comparison.

## Current status (2026-08-08, final)

- AMIP scalar campaign CLOSED; posterior of record =
  `toml/calibration_edmf_fixed.toml` (zsw closure/optics posterior +
  settled microphysics). Neither the EDMF run's displacements
  (under-constrained) nor the panom run's (collapsed, eps_rel +2.0)
  are adopted.
- Next rung: SCM calibrations of the EDMF family (user-side, task
  in progress); fall-speed SCM validation pending.
- Related docs: ROADMAP.md (phases), CALIBRATION_GUIDE.md section 8
  (the ladder), CALIBRATION_LESSONS.md (B9, B10 added from panom).
