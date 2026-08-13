# Session handoff — 2026-08-12 (evening, MDT)

Context for resuming work on branch `kp/ne/amip-calibration` (AMIP
calibration machinery). Written at the end of a local Derecho session;
the next session may be a cloud session (see "Cloud-session limits").

## TL;DR

1. The branch was rebased onto `origin/main` today and ported to
   registered **ClimaCalibrate v0.4.0** (no more branch pin). Validated:
   observations produced by the ported pipeline are **identical** to both
   the pre-rebase code and the original campaign run's vector.
2. The **reproducibility calibration COMPLETED** (2026-08-12 22:15 MDT,
   3 iterations, driver exited 0, workers released). Verdict recorded in
   the repro config header: machinery reproduces exactly (obs and
   iteration-1 parameters bit-identical to the original run), but the
   trajectory diverged because main changed the model physics
   (`cloud_ice_formation: TemperatureDependent`, commit 674ccf36, plus
   ClimaAtmos 0.42.3→0.42.4): G changed 23-41% rel RMS at identical
   parameters. Campaign posteriors must be treated as priors to
   re-verify on the new physics. Comparison tool:
   `compare_trajectory.jl` in the run dir.
3. **Nothing is pushed.** The rebase rewrote history and several key
   files are untracked. A cloud session sees none of this unless it is
   committed and pushed first (`git push --force-with-lease`).

## What is running on Derecho (local-only resources)

- tmux session **`cal_repro`** on the login node runs the calibration
  driver. Reattach: `tmux attach -t cal_repro` (detach `Ctrl-b d`).
  Driver script + env: `/glade/derecho/scratch/kphan/amip_calibration_release_repro_out/driver_repro.sh`.
- PBS worker jobs **7104267, 7104270, 7104271, 7104272** (queue `gpu`,
  account UCIT0011, 12 h walltime, 4 nodes, 13 workers packed 4/node,
  all connected over the HSN at launch).
- Config: `experiments/calibration/amip/config/lwp_clt_swcre_release_repro.jl`
  — byte-identical science to `lwp_clt_swcre_release.jl` (the last
  completed campaign run), but `n_iterations = 3` (half of 6) and output
  to `/glade/derecho/scratch/kphan/amip_calibration_release_repro_out`
  (symlinked from repo root as `amip_calibration_release_repro`).
- Timeline: preflight estimated 1.2 h/iteration, ~4.6 h total → done
  roughly 23:00 MDT 2026-08-12. Worker walltime (12 h) is ample.
- Monitoring rules (README): judge by process liveness and member
  file-write recency; Julia buffers driver.log for hours. `qstat` is a
  cache and can lie for ~30 min.
- If the driver must be killed: check for orphaned `julia-*` jobs in
  `qstat -u kphan` and `qdel` them (guide pitfall 7). The run is
  resumable (rerun the driver with the same CALIBRATION_CONFIG), but
  clean partial member dirs first (preflight checks this).
- The background watchers armed by the previous session died with it.
  Nothing else monitors the run automatically.

## The reproducibility experiment (why this run exists)

Goal: confirm the post-rebase pipeline reproduces the original
floor-release run (`lwp_clt_swcre_release.jl`, run 2026-08-06/07 by
nefrathenrici) before trusting the migrated machinery.

Predictions (also in the repro config header), status:

1. **Observation vector identical to the original run's — CONFIRMED**
   before launch (every sample, covariance block, coverage mask; the
   three-way check original == pre-rebase == post-rebase all matched,
   differences ≤ ~1e-13 relative, i.e. summation-order noise).
2. **Iteration-1 member parameters byte-identical — CONFIRMED, 13/13**
   (`diff` of every `iteration_001/member_*/parameters.toml` against the
   original run's; same priors + seed 42 + TransformUnscented sigma
   points are deterministic).
3. **3-iteration trajectory tracks the original within the weather
   floor — PENDING.** When the run finishes, compare per-iteration
   residuals/loss/parameter means against the original's first three
   iterations. Same-sample comparison is valid: iterations 1–3 grade
   Sep 2006, 2010, 2008 in both runs.

How to judge prediction 3:

- This run's data: `driver.log` steering blocks and
  `iteration_00N/eki_file.jld2` in the repro run dir.
- Original run's data (copied to kphan scratch, member output excluded):
  `/glade/derecho/scratch/kphan/nefrathe_calibration_runs/amip_calibration_release/`
  — has `driver.log`, per-iteration `eki_file.jld2`/`G_ensemble.jld2`,
  plots, and `calibration_report.md`.
- Reference anchors from the original: iteration-1 leverage ratios
  lwp 2.1 / clt 2.6 / swcre 2.4; residual trajectory over its 6
  iterations lwp 0.65→0.66, clt 0.76→0.93, swcre 0.69→0.55 (its σ
  units). Extract per-iteration values from its eki files or driver.log.
- Expected agreement: displacement signs match; residuals within
  ~0.1 σ per iteration. Known confound if it drifts more:
  **ClimaAtmos 0.42.3 (original) → 0.42.4 (now)** — the campaign once
  measured 27% ocean-lwp shift from an Atmos upgrade at identical
  parameters, so a trajectory mismatch indicts the model bump, not the
  calibration machinery (predictions 1–2 already isolate this).
- Also check the iteration-1 go/no-go gate (guide §5): ensemble spread
  of G must beat the noise floor; steering `obs RMS` etc. in driver.log.
- End-of-run: `calibration_report.jl` as a PBS analysis job
  (develop queue, `select=1:ncpus=4:ngpus=1:mem=40GB`), then record the
  verdict in the repro config header (campaign practice), quoting
  numbers, run dir, and dates.

## Repo state (NOT pushed)

- Worktree: `/glade/u/home/kphan/worktree/ClimaCoupler/ne/amip-calibration`,
  branch `kp/ne/amip-calibration`, HEAD `fd733fd0` ("Add AMIP calibration
  machinery") rebased onto `origin/main` = `765727cb`.
- Untracked (commit if they should survive/be visible in cloud):
  - `experiments/calibration/amip/SUMMARY.md` — campaign narrative,
    timeline, loss≠RMSE explainer, obs/covariance changes, epilogue.
  - `experiments/calibration/amip/campaign_reference/` — campaign-branch
    docs (CALIBRATION_LESSONS.md, identifiability_map.md, ROADMAP.md,
    plan docs), analysis scripts (leverage.jl, run_parameter_sweep.jl,
    clt_regime_attribution.jl, …), and `config_with_verdicts/` (all 42
    configs incl. final zsw/edmf/panom verdicts). Provenance:
    nefrathe's clone, branch `ne/calibrate` @ `bfb85870`.
  - `experiments/calibration/amip/config/lwp_clt_swcre_release_repro.jl`
    — the running experiment's config.
  - `HANDOFF.md` (this file).
- The rebase resolutions, in case anything looks surprising:
  - Manifests + Project.toml: main's exactly (ClimaCalibrate 0.4.0
    registered; `[sources]` pin and unused PNGFiles removed).
  - `generate_observations.jl`: OUR statistical machinery kept
    (SVDplusD + noise groups + correlated floors + masks +
    harmonization), ported to the v0.4 SampleBuilder API. Main's version
    had regressed to ScalarCovariance + normalization — deliberately not
    adopted. Port details: `SVDplusDCovariance(; kwargs)` (no date
    positional), one `build_samples_by_times(vars, COVARIANCE_DATE_RANGES;
    FT = Float64)` collection, targets selected by index
    (`findfirst` into covariance dates). **FT = Float64 matters** — the
    new API defaults to Float32.
  - `run_calibration.jl`: merged. Kept `DefaultScheduler(0.1)` (main
    uses DataMisfitController, which the campaign showed collapses
    ensembles — do not "fix" this), CALIBRATION_CONFIG env selection,
    obs vector in output_dir, worker packing 4/node, HSN-bind exename
    hook, 11 h empty-pool timeout, priors plot. Adopted main's
    TEST_CALIBRATION-on-CPU and `mem = "16GB"` for the test path.

## Environment

- Modules/Julia: `module load climacommon/2025_02_25`;
  `JULIA=/glade/campaign/univ/ucit0011/software/julia/julia-1.11.3/bin/julia`;
  project `--project=experiments/AMIP`.
- Versions: ClimaCalibrate 0.4.0 (registered), EKP 2.7.1, ClimaAnalysis
  0.5.23, ClimaAtmos 0.42.4. Freeze while the run is in flight.
- Login nodes: 10 GiB cgroup → `--heap-size-hint=3G` and
  `JULIA_NUM_PRECOMPILE_TASKS=8` (parallel precompile workers otherwise
  get OOM-killed, which once left EMPTY package dirs in
  `~/.julia/packages/<Pkg>/<slug>` — instantiate then thinks the package
  is installed and precompile fails with "does not seem to be
  installed"; fix is deleting the empty slug dirs and re-instantiating).

## Cloud-session limits

A cloud session cannot: reach Derecho scratch/tmux/PBS, monitor or
manage the run, or see unpushed/untracked work. It CAN: work on the
committed+pushed repo (docs, code review, analysis scripts, planning).
Run monitoring, the trajectory comparison, the report job, and any
future launches need a session with Derecho shell access.

## Other completed work this session (context)

- Fixed the broken environment that predated the rebase: Manifest
  pinned an unfetchable tree of `ClimaCalibrate#ne/worker-packing`
  (never-pushed or force-pushed commit) → resolved to branch head; then
  the empty-husk depot problem (see Environment above).
- Session artifacts (scratchpad, ephemeral — may not survive):
  observation A/B snapshots and compare scripts under
  `/glade/derecho/scratch/kphan/tmp/claude-48097/.../scratchpad/`
  (`snapshot.jl`, `compare_snapshots.jl`, `obs_pre/`, `obs_post/`,
  `snapshot_original_release.jld2`). The methodology: dump samples,
  cov diagonals + Frobenius/sums, masks to base-type JLD2, compare
  across environments.
- Campaign endgame (recovered from ne/calibrate, now in SUMMARY.md
  epilogue and campaign_reference/): zsw (2026-08-07) halved the
  clt-for-swcre trade; edmf (2026-08-08) flat — AMIP-invisible, scalar
  campaign converged; panom (2026-08-08) pattern residual immobile,
  ensemble collapsed — clt pattern error is directionally orthogonal to
  the entire parameter space; campaign closed toward model development /
  the SCM rung. Attribution: residual concentrates in the Sc-to-Cu
  transition (clt_regime_attribution.jl).
