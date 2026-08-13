# Hard-won calibration rules

Rules that this campaign learned from runs that gave the wrong answer, or
that looked correct and were not. Each rule carries the evidence that
produced it, because a rule without its evidence gets argued away later.

CALIBRATION_GUIDE.md describes how to run a calibration. This file
describes how to avoid being misled by one. The two overlap in places and
should be merged when someone has time.

Status: collected 2026-08-04, not yet merged into the guide. Two proposed
preflight and steering changes are listed at the end and are NOT
implemented.

## A. Reading the numbers

### A1. Spread contraction is not identifiability

Use mean displacement in prior-sigma units. A parameter that ends at its
prior mean learned nothing, however tight its ensemble looks.

Evidence: every one of vterm's 6 parameters contracted to 0.042-0.068x
prior spread, and every one of relax's 8 to 0.03-0.04x, whether or not the
data touched them. That uniformity is the fixed-timestep scheduler's
accumulated pseudo-time, not information. In the same vterm run, mean
displacement separated the parameters cleanly: v_snow +1.45 sigma,
rain_tau -0.71, q_liq -0.46, v_rain +0.42, subl_dep +0.39, v_ice -0.06.
Only v_ice was untouched, and only mean displacement showed it.

Consequence for the record: any earlier verdict that cited a spread ratio
as proof of identifiability needs rechecking. relax's prediction 2 about
cloud_fraction_floor_release_abs_margin is one of them: it cited a
0.03-0.04x spread, but the margin's mean moved 1.0 to 1.066, about 0.13
sigma, which is barely at all.

### A2. Check an observable's whitened misfit before trusting it

An observable whose misfit at the prior is well below 1 contributes no
gradient. The calibration will ignore it, and no amount of iteration will
change that.

Evidence: in the clt-only run, clt's whitened misfit at iteration 1 was
0.122, so its RMS residual was 0.35 of its assumed noise. Over two
iterations clt went 0.620 to 0.623, that is, slightly worse, with nothing
competing for the parameters. The optimizer was correct to do nothing.

How to compute it: at iteration 1 the prior misfit term is exactly zero,
so `EKP.get_error(ekp)[1] * n_params` is the whitened misfit
`(1/d) (y - g)' inv(Gamma) (y - g)`. Its square root is the RMS residual
in units of the assumed noise.

### A3. Diagonal RMS and full-covariance misfit can disagree completely

Report both. The diagonal number describes the size of the error. The
whitened number describes whether the error has a shape the noise model
expects.

Evidence: between vterm iterations 3 and 5 the point-count-weighted
diagonal mean square residual was flat, 0.420 to 0.418, while the full
misfit rose 29 percent, 0.181 to 0.234. The residual amplitude did not
change. Its pattern moved into directions that the correlated covariance
considers improbable. Reading only the diagonal RMS called a degrading fit
stable.

### A4. Know exactly which error metric you are looking at

`EKP.get_error` returns `bayes_loss` when `impose_prior = true`, which is
`misfit / n_params + prior_misfit / n_obs`. `EKP.Visualize` with
`error_metric = "loss"` returns the data misfit. They differ by a factor
of the parameter count, and the prior term is negligible when n_obs is
large. Confirm the relation by arithmetic before comparing runs: 0.0302 x
6 = 0.181 matched the plotted value exactly.

### A5. Quote the posterior at the loss minimum when the loss is not monotone

Say which iteration you used. vterm's misfit fell to a minimum at
iteration 3 and then rose to its starting value by iteration 5. relax
shows the same shape on the same observations, so it is a property of the
fixed timestep against an over-informative observation, not of the
parameter set.

## B. Designing the experiment

### B1. The correlated noise floor decides which errors you may care about

Choose the decorrelation length shorter than the error structure you want
to penalize.

Evidence: with an 800 km floor, clt's whitened residual was 0.35 of the
assumed noise while its diagonal residual was 0.62. The gap is the
correlated floor absorbing the smooth part of the error. A near-global
cloud-cover deficit is exactly the kind of error such a floor makes free,
and it is also exactly the error we wanted to fix.

### B1a. Leverage is what limits a calibration, and no weighting buys it

Measure `residual / ensemble spread` per observable at iteration 1 (see
`leverage.jl`). It is the excursion, in units of the current ensemble
width, needed to close the residual. This ratio is **invariant to the noise
model**: both terms are divided by sigma, so scaling `model_error_scale`
cancels exactly. A permissive floor makes an observable look already fit; a
tight floor makes it look badly fit; neither changes whether the parameters
can reach it.

Consequence: when the ratio is high, the fix is a parameter the observable
responds to. It is never a reweighting, and it is not a smaller region.

Evidence: clt's ratio was 3.7 globally with the cloud-fraction shape frozen
and 2.7 with it free, and clt improved only in the second case. Narrowing to
a 45S-65S band to concentrate the signal made every ratio WORSE, lwp 1.5 ->
3.7, clt 3.7 -> 5.7, swcre 2.6 -> 3.6, and three iterations moved the
residuals about 2 percent. Measured anchor points for the ratio are in the
`leverage.jl` header.

### B1b. Never launch a regional run before subsetting an existing one

Any question of the form "would this work better in region X" is answerable
for free from a G ensemble you already have. Restrict an existing run's
observation vector to the latitudes of interest and recompute
residual/spread there. Only launch the regional run if that number is
promising.

Evidence: the global vterm run at iteration 1, restricted to 45S-65S, gives
ratios of 3.7 (lwp), 5.7 (clt) and 3.7 (swcre). The dedicated Southern Ocean
band run, which differed in parameter set (q_liq fixed), grid (2.5 instead
of 5 degrees) and floor (400 instead of 800 km), measured 3.7, 5.7 and 3.6.
Identical. Those ~35 GPU-hours bought a number already on disk.

The same subset also tells you which region WOULD work. For this parameter
set the answer is the tropics, not the Southern Ocean:

    band       lwp   clt   swcre
    global     1.5   3.7   2.6
    45S-65S    3.7   5.7   3.7
    30S-30N    1.4   3.4   2.6
    30N-65N    1.5   3.7   2.0

Read the spread column, not the residual, to see why. Parameter response in
the Southern Ocean is 3-4x weaker than in the tropics (lwp 0.173 against
0.582) while the biases are comparable or smaller. swcre is the clearest
case: the model is LEAST biased there (0.360 sigma against 1.077 in the
tropics) and simply does not respond. These are warm-rain and fall-speed
parameters, and they act where there is thick liquid cloud and active rain.

### B1c. Spatial reduction buys signal-to-noise, and the raw ratio hides it

Collapsing longitude to a zonal mean removes incoherent weather while
keeping the coherent parameter response. Measured on vterm's own G
ensembles, in physical units, with the weather floor taken from the
iteration-6 spread and the signal by quadrature subtraction:

    obs          weather floor   parameter signal    S/N
    lwp 2-D             0.0168             0.0231   1.37
    lwp zonal           0.0049             0.0116   2.38
    clt 2-D             0.0437             0.0318   0.73
    clt zonal           0.0121             0.0241   1.99
    swcre 2-D             6.11               5.85   0.96
    swcre zonal           1.25               3.53   2.83

Reduction improves S/N by 1.7 to 2.9x and lifts clt and swcre above 1 for
the first time. Reachability is essentially unchanged: noise-corrected
ratio goes lwp 1.9 -> 2.6, clt 7.4 -> 7.5, swcre 4.1 -> 4.4.

**The trap:** plain residual/spread gets WORSE under averaging (lwp 1.6 ->
2.4) because averaging removes weather from the denominator. Whenever the
reduction changes, divide by the parameter signal, not the total spread, or
you will conclude that cleaner data is worse data. This is a limitation of
the ratio reported by `leverage.jl`.

Two things break when the spatial dimension collapses, both worth knowing
before trying it:

- `QuantileRegularization(qtl)` needs at least `1/qtl` values per variable.
  The 0.05 default needs 20; a zonal mean over a 20 degree band leaves 8
  and errors with "Insufficient samples for computing quantile". Set
  `OBS_REGULARIZATION_QUANTILE` to about `2/n_values`.
- A floor set by `model_error_scale` still scales with field magnitude,
  while averaging shrinks the field's spatial structure. swcre's zonal
  noise came out at 132 percent of its signal RMS at the same 0.75 that
  gave 71 percent on the 2-D grid. Re-tune floors after reducing.

### B2. Ask for the minimum run that separates the hypotheses

It is often one run, not two. icabl_sub answered a two-arm question by
matching its data to an existing run, so the second arm was never needed.

### B3. Freezing parameters at a posterior is cheap and it localizes function

Cost falls with the parameter count, since the ensemble is 2p+1. Freezing
also tells you what each frozen parameter was doing.

Evidence: vterm froze eps_rel, steepness, the floor-release margin,
snow_tau and cond_evap at relax's posterior and opened the three fall
speeds. Only clt degraded (0.620 to 0.654 across iterations) while lwp and
swcre improved. That identified the cloud-fraction shape parameters as
clt's levers more directly than a sensitivity study would.

### B4. Put the decisive number in the first update, not the last

Then cancel when the answer arrives. The clt-only run was answerable at
iteration 2 and was cancelled at 2 of 6, saving about 65 GPU-hours.

Related: identical prior, seed and ensemble construction give a
bit-identical iteration 1 across runs, because TransformUnscented's sigma
points depend only on the prior and the dimension and the forward model is
deterministic. Iteration 1 is therefore a free cross-run consistency check
and never carries new information.

### B5. Write falsifiable predictions before launch, with a metric and a threshold

Be willing to conclude that the metric was wrong. That is a real outcome
and often the most valuable one.

Evidence: rule A1 exists only because vterm's prediction 2 was falsified
on its stated metric (ice spread stays above 0.5x prior). The intended
claim, that these observables carry no ice information, was correct. The
metric was not. A conclusion written after the run would have recorded
neither.

**B9. Leverage ratios built from response magnitudes over-promise;
direction decides.** The pattern-anomaly scoping (2026-08-08) computed
ratio = residual / quadrature-sum-of-responses ~ 1 for clt's zonal
anomaly, and the calibration built on it (panom) moved the residual not
at all (1.41 -> 1.40 sigma, same-Sep endpoints) while eps_rel ran +2.0
prior sigma and the ensemble collapsed to 0.0x spread. Parameters moved
the pattern strongly, but not TOWARD the observed pattern: magnitude
without alignment. Before buying a run from a scoping ratio, project
the per-parameter responses onto the residual vector and gate on that
projection. (The production iteration-1 leverage in leverage.jl shares
this limitation: its "spread" is magnitude, not alignment.)

**B10. Anomaly-type observables over-inform under the per-point floor.**
The floor is (scale x seasonal_mean_i)^2 per point; an anomaly's means
pass through zero, so large regions get a near-zero floor and the block
carries far more nominal information than the field justifies. panom
collapsed by mid-run under the same mechanism as the 2-D diagonal-floor
runs (B-series). A metric that transforms a field toward zero mean needs
its own floor design (for example scale x field RMS as a uniform floor),
not the default.

## C. Operations

### C1. Run directories must be symlinks into scratch

A real directory in the repo root writes into a 100 GiB home quota that
normally sits near 90 percent. The failure is not a clean "disk full":
JLD2 checkpoint writes truncate, every member throws EOFError from
`Checkpointer.checkpoint_model_cache`, and the driver dies. Create the
scratch directory and the symlink in the same step, and confirm `ls -ld`
shows an arrow before qsub.

### C2. The driver log's final `exited` line is ground truth, not qstat

PBS query outages reported three separate false vanished-job alarms during
this campaign while the drivers were running normally. Use qstat only to
detect a duplicate driver on one run directory.

### C2a. Two jobs sharing one `PBS -o` path clobber each other

If you fix a script and resubmit while the first job is still queued, the
first job runs the copy PBS took at submission time and writes its old
failure to the same log. Reading that log then shows a bug you already
fixed. Wait on the QUEUE being empty of that job name, not on the log
having content. Cost this session: two wrong diagnoses in ten minutes.

### C3. `PBS -o` appends across resubmissions

Scope every log grep to the current run, for example with an awk range
anchored on the submit script's start banner. An unscoped grep diagnosed a
crash from the previous attempt and reported it as a live failure.

### C4. Detect dead workers by cput growth, not by job state

Compare `resources_used.cput` growth against elapsed time per node.
relax_amip lost four members to three dead workers on one node while PBS
reported all four node jobs as running; that node showed 4:13 of cput over
4:13 of walltime and 22 GB of memory, against about 16:50 and 78 GB on its
healthy peers.

## Open items, not implemented

1. Preflight check: compute each observable's whitened misfit at the prior
   and FAIL below about 0.5. This would have caught clt before two runs
   were spent on it. See A2 for the one-line formula.
2. Steering indicator: report mean displacement in prior-sigma units
   beside the spread ratio, and flag any parameter that ends within 0.1
   sigma of its prior mean. See A1.
3. Merge the overlapping parts of this file into CALIBRATION_GUIDE.md.
4. Recheck every verdict that used a spread ratio as identifiability
   evidence, starting with relax's floor-release margin (A1).
