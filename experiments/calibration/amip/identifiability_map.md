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

**Model version: ClimaAtmos 0.42.2** (all numbers below; older runs on 0.41.3 are
not comparable).

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

## Rows 3+ — planned (ROADMAP.md priority order)

3. **Cloud optics / droplet number / effective radius** — targets swcre's "too
   bright" axis at ~fixed lwp.
4. **Cloud fraction scheme** (quadrature/eps params) — targets cl directly.
5. **Ice microphysics** — targets lwcre/clivi. **Mixing length / diffusion** —
   targets ta/hur pressure-level profiles.

---

## Current status (2026-07-28, 04:10 MDT)

- **Row-2 sweep**: ✅ COMPLETE (jobs 6930400–6930412, 13/13 members). Results and
  verdicts in the Row-2 section above. Two follow-ups queued for a future
  mini-row: `mixing_length_static_stab_coeff` / `mixing_length_Ri_crit` (the wired
  TKE-dissipation controls), and a combined-detrainment point (dmfvd↑ + detr_buoy
  at base) to test whether the family can close more of the swcre gap jointly.
- **lwp+pr calibration** (job 6926553): iteration 8 of 10, q_liq = 2.80e-4,
  rain_τ = 1512 — stable; walltime ends ~06:20 MDT, expect ~iteration 10.
- **SON spin-up diagnosis**: CONCLUDED (see Replication note in Row 1) and stopped
  at 3/9 iterations — the continuation manager hit a **restart-diagnostics bug**:
  members resumed mid-iteration emit wrong-dated monthly averages (diagnostics
  windows restart at the model checkpoint), observation_map crashes with "Too many
  NaNs", and every resubmission re-enters the wedge. Both arms had already
  converged, so no further GPU is spent. Bug recorded; PACKAGING_PLAN resume
  semantics must clean partially-run members (see plan addendum).
- Related docs: ROADMAP.md (phases), PACKAGING_PLAN.md (freeze framework),
  STEERING_INDICATORS_PLAN.md (log-only indicators, v1 scoped).
