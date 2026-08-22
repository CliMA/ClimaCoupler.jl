# Calibration handoff — 2026-08-21 (night, MDT)

How to run and babysit AMIP-mode calibrations on this branch
(`kp/ne/amip-calibration`), written for an incoming agent. Supersedes the
2026-08-12 repro-run handoff (see git history for it). Campaign results and
post-mortems live in the config headers and the session memory; this file is
the OPERATIONAL knowledge.

## 0. THE JOB AT HAND: rlut_son_ocean (LIVE — read this section fully)

**Identity**: AMIP-mode, ocean-only, SON-2010 seasonal-mean rlut calibration.
Config `experiments/calibration/amip/config/rlut_son_ocean.jl` (its header
holds the rationale and four pre-registered predictions — grade them at the
end). 7 members, target **4 iterations** (T=0.4), n = 1696, floor 4.5 W/m²,
expected iteration-1 whitened residual ≈ 1σ.

**Run dir** (everything lives here):
`/glade/derecho/scratch/kphan/amip_calibration_rlut_son_ocean_out`
(= repo symlink `amip_calibration_rlut_son_ocean`). Members are ~122
simulated days ≈ **4.3 h wall each**; an iteration = one member's wall time
(7 members run concurrently on 2 GPU nodes).

**What is running**: tmux session `cal_son_chain` executing
`<run_dir>/chain_son.sh`, started 2026-08-21 21:10 MDT:

    step 1  wait for the Aug-1 atmosphere IC            DONE (21:10)
    step 2  one-week smoke test                          verdict ~22:00 Aug 21
    step 3  launch 1: CALIBRATION_N_ITERATIONS=2         done ~07:15 Aug 22
    step 4  launch 2: CALIBRATION_N_ITERATIONS=4         done ~16:00 Aug 22

Any failure STOPS the chain with a diagnosis line. It never retries by
itself — that is deliberate (see the restart hazard, §3).

### Status check (the only command you need routinely)

    cat /glade/derecho/scratch/kphan/amip_calibration_rlut_son_ocean_out/chain.log
    # supporting evidence:
    tmux ls                                  # cal_son_chain present?
    qstat -u kphan                           # smoke job, or 2x julia-* workers
    find <run_dir> -maxdepth 2 -name G_ensemble.jld2 | wc -l   # iterations done

### Scenario playbook

**A. Chain healthy** → touch nothing. Progress markers per stage: smoke →
`SMOKE TEST PASSED` in `/glade/derecho/scratch/kphan/rlut_son_smoke/smoke.log`;
launches → `Running member N` lines in `<run_dir>/driver.log`, member
`output.log`s updating (they print progress + ETA), `iteration_00N/` dirs
appearing with `G_ensemble.jld2` when an iteration closes.

**B. Chain stopped: "smoke test FAILED"** → read the smoke log's ERROR +
stacktrace. Do NOT launch the calibration. Likely classes: IC file/date
resolution (see §5 AMIP-IC fact), AMIP-mode config validation, instability
(look for `Found NaN` / `Y.f.sgsʲs`). Fix, rerun the smoke by hand:

    cd <repo>; mkdir -p /glade/derecho/scratch/kphan/rlut_son_smoke
    qsub -v REPO=$PWD,CALIBRATION_CONFIG=$PWD/experiments/calibration/amip/config/rlut_son_ocean.jl,SMOKE_DIR=/glade/derecho/scratch/kphan/rlut_son_smoke,SMOKE_SIM_SECONDS=604800 \
         -o /glade/derecho/scratch/kphan/rlut_son_smoke/smoke.log \
         experiments/calibration/amip/smoke_test.sh

**C. Chain/tmux died mid-launch (login node reboot, kill)** — the one case
needing care:
  1. `qdel` any surviving `julia-*` worker jobs FIRST (never run two drivers
     against one pool).
  2. Count completed iterations (`G_ensemble.jld2` count = K).
  3. **If an iteration is INCOMPLETE** (an `iteration_00(K+1)` exists with
     member dirs but no `G_ensemble.jld2`): delete its member subdirectories
     before relaunching —
     `rm -rf <run_dir>/iteration_00(K+1)/member_*` — so members rerun from
     scratch. Otherwise `detect_restart_files` resumes them from mid-month
     checkpoints and their monthly means are silently wrong (§3).
  4. Relaunch cleanly, target sized to the 12 h worker walltime
     (2 iterations per launch MAX for this run):

         tmux new-session -d -s cal_son_l2 \
           'CALIBRATION_N_ITERATIONS=<K+2 (cap 4)> bash /glade/derecho/scratch/kphan/amip_calibration_rlut_son_ocean_out/driver_rlut_son_ocean.sh'

     Resume skips the K completed iterations automatically. Repeat until 4.

**D. A launch exited but below its target** (chain says STOP with N/target)
→ `grep -E "Error|Found NaN|failed" <run_dir>/driver.log | tail`. Member
NaNs: EKP absorbs ≤3/7 dead members but note WHO died (σ=0.3 priors were
audited stable — a death is news; record it). 100% failure aborts the driver
("Execution halted") — that is a config/physics problem, not a retry case.

**E. All 4 iterations done** → post-run sequence:
  1. Gate + spread:
     `GNG_MIN_SPREAD=1.5 julia --startup-file=no --project=experiments/AMIP experiments/calibration/amip/go_no_go.jl <run_dir> 1`
  2. **Physical trajectory** (the gate misses this — lwcre lesson): per
     iteration, load `iteration_00N/G_ensemble.jld2`, drop all-NaN member
     columns, mean over members, compare to
     `EKP.get_obs(JLD2.load_object("<run_dir>/observation_vec.jld2")[1])`:
     bias and RMS per iteration on identical weather. Bias should shrink;
     RMS flat is expected; BOTH degrading while EKP loss falls = the
     pattern-steering pathology, stop trusting the posterior direction.
  3. Parameters: `iteration_005/eki_file.jld2` + `iteration_001/prior.jld2`,
     `EKP.get_ϕ(prior, ekp)` per iteration. Key questions, from the config
     header's predictions: all 7 alive in iteration 1? whitened residual
     ≈1σ? contraction <30x? and THE question — does c6 rise again, and if
     so does lwcre (free in the saved rlutcs/rlut diagnostics of every
     member) degrade in mirror (compensating error) or hold (real signal)?
  4. Plots auto-generate per iteration (`bias_sample_dates.png`,
     `g_vs_obs.png`, first/last strips in the run dir). The bias-map panels
     are MEMBER 1 (= the UKI mean member), not the ensemble mean.

### Launch anatomy (what the driver does, for debugging)

`driver_rlut_son_ocean.sh` (in the run dir): loads climacommon, exports
`CALIBRATION_CONFIG`, `CALIBRATION_WORKER_EXENAME` (HSN bind — required for
login-node drivers), `CALIBRATION_N_ITERATIONS`; runs
`run_calibration.jl` under `--heap-size-hint=3G`, tees to `driver.log`.
The driver loads `observation_vec.jld2` (hard-fails if absent — regenerate
with the prep, §2), builds the EKP object, submits 2 packed GPU worker jobs
(12 h walltime, 4 workers/node), farms members, updates, writes
`iteration_00N/` per iteration, exits after `CALIBRATION_N_ITERATIONS`.
Worker PBS stdout persists in `<run_dir>/worker_logs/`.

## 1. Anatomy of a calibration

Selection is by environment variable: `CALIBRATION_CONFIG=<config .jl>`.

| layer | example (current run) | notes |
|---|---|---|
| calibration config | `experiments/calibration/amip/config/rlut_son_ocean.jl` | loss variables, dates, priors, noise model, output_dir, flags (`OCEAN_ONLY`, `SEASONAL_MEAN`); headers carry run rationale + pre-registered predictions — READ THEM |
| coupler YAML | `config/amip_configs/amip_calibration_pigroups_son.yml` | mode, grid, start_date, `coupler_toml`, `edmfx_entr_model`, `log_to_file` |
| atmos YAML | `config/atmos_configs/climaatmos_progedmf_1m.yml` | merged FIRST; coupler YAML keys override it |
| parameter TOML | `toml/amip_progedmf_1m_pigroups.toml` | the atmos YAML's own `toml:` is REPLACED, not layered — only `coupler_toml` + the member file reach ClimaParams |
| member file | `iteration_N/member_M/parameters_spliced.toml` | `<base>_E<i>` prior names splice into vector params (`entr_param_vec`); base vector MUST exist in a coupler_toml |

Environment: `module load climacommon/2025_02_25`, Julia
`/glade/campaign/univ/ucit0011/software/julia/julia-1.11.3/bin/julia`,
`--project=experiments/AMIP` (manifest pins ClimaAtmos `#main`). Ground truth
of what a member ran: its resolved `*_parameters.toml` (341 params) and `.yml`
in the member's `clima_atmos` output dir.

Run dirs live on scratch, symlinked from the repo root
(`amip_calibration_<name>` -> `/glade/derecho/scratch/kphan/..._out`).
The observation vector is keyed to `output_dir` — new run = new dir, or you
silently grade against a stale covariance.

## 2. The procedure (gates in order; skip none for a new configuration)

1. **Prep job** (develop queue, ~8 min; template: `<run_dir>/prep.sh`):
   instantiate → `generate_observations.jl` → `plot_observations.jl`
   (eyeball `observation_check.png` — reconstructs what EKP actually sees)
   → `check_g_dates.jl` (synthetic member through the REAL preprocess +
   GEnsembleBuilder; catches obs/sim date misalignment at zero GPU)
   → `preflight.jl` (wiring check D is the one that matters: a prior the
   model silently ignores must fail here, not appear as a flat posterior).
   PASS = `0 failed`, all `wiring` PASS, and **zero `KeyError` strings in
   the log** (the default logger swallows message-construction errors).
2. **Sigma-point audit**: print the actual constrained values the 7 members
   will run (see §5 stability envelope). One minute on a login node.
3. **Smoke test** (`smoke_test.sh`): `SMOKE_SIM_SECONDS=360` = construction
   +2 steps (~15 min); `604800` = one week; days mode via `SMOKE_DAYS`.
   Mandatory for any new mode/physics; it has caught a fatal config error
   or instability before every launch that needed one.
4. **Launch** in tmux on a login node (never a compute job for the driver).
5. **Iteration-1 go/no-go** (§4).

## 3. Walltime and the restart hazard (the most important operational rule)

Workers are PBS jobs capped at **12 h** (queue max). Members killed at
walltime and rerun are DANGEROUS: `detect_restart_files=true` restarts from
checkpoints, but monthly-diagnostic accumulation state is NOT checkpointed,
and restarts may write to a fresh `output_NNNN` segment while the observation
map reads `output_active` only → silently wrong or missing months in G.

**Rule: never let walltime kill members. Size launches to clean iteration
boundaries** via `CALIBRATION_N_ITERATIONS` (env-overridable in the config)
so each driver run exits cleanly inside 12 h, then relaunch with a higher
target (resume skips finished iterations). Throughput ≈ 30 simulated
days/GPU-hour at h_elem 12; iteration wall time = one member (all members run
concurrently on 2 packed GPU nodes = 8 slots).

Mitigations already in place: member checkpoints are month-aligned
(`model_interface.jl`), and `check_season_months` errors loudly on missing
months. For unattended multi-launch there is `calibration_relay.sh
<run_dir> <target> <driver cmd...>` + auto-started watchdog (tested), and the
chain-script pattern (`chain_son.sh`: wait-for-file → smoke → launches).

## 4. Reading results (and the trap)

- `go_no_go.jl <run_dir> <iter>`: leverage ratio (≤3), physical reachable
  spread (`GNG_MIN_SPREAD`, scale it to the run's floor), contraction (≤10x),
  dead members. Exit 1 = stop the run.
- **The gate is not sufficient** (lwcre lesson): EKP minimizes the
  Σ⁻¹-weighted misfit and can "improve" while the loss variable's physical
  bias/RMS degrade. After iterations 1→2, compute same-weather bias/RMS of
  the loss variable from `G_ensemble.jld2` vs the observation; two
  consecutive wrong-direction moves = stop.
- Loss/params live in `iteration_<last>/eki_file.jld2` (cumulative: full
  history) + `iteration_001/prior.jld2`. `EKP.get_error` across a
  failed-member boundary is not comparable. Contraction reflects the ASSUMED
  noise, not achieved fit — 24.8x–7000x observed under tight floors; treat
  posterior spreads as unreliable, means as the result.
- Logs: `driver.log` (scheduling, EKP, thrown exceptions incl. member NaNs),
  member `output.log` (live, `log_to_file`), `worker_logs/worker-*.out`
  (persistent PBS stdout), `relay.log`/`chain.log`. qstat lags ~30 min on
  login nodes — judge liveness by file mtimes, never one qstat miss.

## 5. Hard-won facts (violate any of these and you rediscover a post-mortem)

- **model_error_scale is per-field**: floor = (scale × field mean)²; the
  honest criterion is floor ≈ interannual spread. Measured: rlut ocean 0.02;
  lwcre honest ≈ 0.2 (0.02 gave 7σ whitening and one-update collapse). Never
  transplant a scale between fields. Too-small floors also steer updates
  into near-null pattern directions (optimizer diverges from physical
  metrics).
- **Prior stability envelope** (PiGroups): c2/c3 sigma points at ±0.87
  NaN'd (`Y.f.sgsʲs` updraft blow-up); ±0.52 (σ=0.3) is inside. Losing
  members in iteration 1 corrupts the DECISIVE update via failure-handling
  artifacts. Audit sigma points before every launch.
- **SEASONAL_MEAN**: uses `ClimaAnalysis.average_season_across_time`
  symmetrically; requires MONTH-ALIGNED start dates (a partial leading slab
  would contaminate a season silently); `check_season_months` guards the
  requested windows on both sides.
- **AMIP vs subseasonal ICs**: `initial_condition: WeatherModel` is the
  coupler DEFAULT even in amip mode; the atmos IC file
  `era5_init_processed_internal_<start>_0000.nc` resolves from
  `era5_initial_condition_dir` (atmos-only in amip mode) else the
  `wxquest_initial_conditions` artifact, which has FIXED dates — check
  before choosing a start_date. Subseasonal mode needs all seven processed
  IC products; amip needs only the one atmos file.
- **All coupler TOMLs need `type = "float"` fields** — ClimaParams' override
  warning indexes `entry["type"]` and a missing field KeyErrors every member
  (still unfixed upstream).
- Login nodes: `--heap-size-hint=3G` (10 GiB cgroup), `--startup-file=no`
  (user startup.jl is broken under this project), `JULIA_NUM_PRECOMPILE_TASKS`
  ≤2 when instantiating on login. Artifacts resolve via
  `~/.julia/artifacts/Overrides.toml` → symlink to
  `ClimaArtifacts2/artifacts/Overrides.toml` (if CERES/radiation_obs "not
  found", the symlink is missing).
- Workers: submitted by the driver (`julia-*` job names), inherit env via
  `qsub -V`, must bind the HSN via `CALIBRATION_WORKER_EXENAME=
  experiments/calibration/amip/julia_worker_hsn_bind.sh` when the driver is
  on a login node. `DefaultScheduler(0.1)`: T = 0.1 × iterations; 4
  iterations = a deliberate 40% posterior.

## 6. Campaign state (details in config headers + session memory)

- micro_edmf (8 params, lwp/swcre/lwcre): null; EDMF params flat.
- rlut_pigroups (land+ocean, scale 0.05): null — floor = total error
  declared the residual irreducible; ≤4% of it was reachable anyway.
- rlut_pigroups_ocean (scale 0.02): bias 93% removed, c6 0.30→0.34, BUT
  shown to be a compensating-error fit (clear-sky Δrlutcs ≈ −6.5 W/m²;
  lwcre degraded monotonically as c6 rose).
- lwcre_pigroups_ocean (scale 0.02 = mispriced): one-update collapse
  (7000x), c6 up again via pattern-steering while lwcre's own bias/RMS
  degraded. Killed at iteration 4 by choice.
- **rlut_son_ocean (running)**: AMIP mode, SON-2010 seasonal mean, ocean,
  σ=0.3 priors, scale 0.02 (floor/interannual ≈ 1.2 — in band). Predictions
  pre-registered in the config header; grade them.
- Queued ideas: JOINT rlut+lwcre (+rlutcs) loss with per-variable scales;
  one-at-a-time perturbation gate before any new parameter set; upstream
  PRs (ClimaParams `get(entry,"type",…)`; type fields on main's TOML;
  worker exceptions logged member-side; singleton-dim squeeze in plotters).
