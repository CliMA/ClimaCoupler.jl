# Packaged Calibration Framework (freeze-from-snapshot)

> STATUS (2026-07-28): ON HOLD by user decision. Phase-1 snapshot/provenance/
> fingerprint/verify functions exist in `experiments/calibration/amip/packaging.jl`
> (uncommitted, untested). The launcher and driver template are not built.

## Context

A calibration's definition is currently scattered across five layers: the config `.jl`
(globals: PRIORS, COVARIANCE_DATE_RANGES, levels…), the driver scripts (hard-coded EKI
hyperparameters: `DefaultScheduler(0.1)`, `model_error_scale`, regularization), the
coupler/atmos YAML + param TOML (referenced by mutable repo-relative path), the Julia
environment (`experiments/AMIP/Manifest`), and the shell (ENV["CALIBRATION_CONFIG"],
hand-written PBS scripts loose in scratch). **Nothing environmental is recorded or
enforced** (verified: ClimaCalibrate has zero provenance capture).

Motivating incident: a mid-study `git rebase origin/main` moved experiments/AMIP from
ClimaAtmos 0.41.3 → 0.42.2, changing simulated ocean LWP by **27% at identical
parameters** — silently invalidating five runs' worth of cross-run comparison
(`q_liq≈7e-4` was obsoleted by it). This must be structurally impossible.

User decisions: **freeze semantics** (run executes from an immutable snapshot; later
repo edits can't touch a running/resumed run), **full lifecycle** (obs gen → calibrate →
resume → report, + parameter sweep), **lives in ClimaCoupler** (experiments/ +
src/CalibrationTools.jl).

## Key facts the design exploits (verified)

- `experiments/AMIP/Manifest-v1.11.toml` resolves ClimaCoupler via **`path = "../.."`**
  → copying the repo tree makes the copied env resolve ClimaCoupler to the copy
  automatically. No Manifest surgery.
- All other deps (incl. ClimaCalibrate `ne/async`) are pinned by git-tree-sha1 to
  **immutable, content-addressed** dirs in `~/.julia/packages` → a copied Manifest
  reproduces the exact dependency code without network.
- The git-tracked tree is only **18 MB** → full-repo snapshot per run is trivial.
- Scripts locate everything via `pkgdir(ClimaCoupler)` → when run under the snapshot's
  project, `pkgdir` points INTO the snapshot, so config/toml/YAML paths freeze for free.
  **Exception**: config files define `output_dir = joinpath(pkgdir(ClimaCoupler), ...)`,
  which would put output inside the snapshot — output location must come from the
  launcher instead.
- ClimaCalibrate restart = `last_completed_iteration(output_dir)`; resuming is just
  re-running the driver with the same output_dir.
- ClimaCalibrate's WorkerBackend workers inherit env via `qsub -V` +
  `JULIA_PROJECT` propagation from the active project → workers automatically use the
  snapshot project if the driver does.

## Run-package layout

```
/glade/derecho/scratch/nefrathe/calibration_runs/<name>_<YYYYmmdd_HHMM>/   (= output_dir)
├── package/                          # immutable after launch
│   ├── code/                         # git archive of ClimaCoupler @ SHA, dirty patch applied
│   │   ├── src/  config/  toml/  experiments/ ...
│   ├── provenance.toml               # see schema below
│   ├── dirty.patch                   # uncommitted changes at launch (may be empty)
│   └── launch.pbs                    # generated driver script (self-contained)
├── observation_vec.jld2              # generated FROM the snapshot at launch
├── coverage_masks.jld2
├── iteration_NNN/…                   # ClimaCalibrate layout (eki_file, G_ensemble, members)
├── calibration_report.md
└── driver.log, worker_*.log
```

`output_dir` root stays exactly where ClimaCalibrate + existing scripts expect their
files (observation_vec at root, iteration_NNN under it) — no layout API changes.

## provenance.toml schema

```toml
[launch]   time, host, user, julia_version, julia_binary, run_dir, mode  # calibrate|sweep
[git]      sha, branch, describe, dirty(bool), dirty_patch_sha256
[environment]
  manifest_sha256, project_sha256
  # extracted for at-a-glance drift detection (the Atmos incident would be one line):
  climaatmos_version, climacore_version, climacalibrate_rev, climacoupler_version
[spec]     config_path, config_sha256, coupler_yaml_sha256, atmos_yaml_sha256, toml_sha256
[observations]  observation_vec_sha256, coverage_masks_sha256, generated_at
[external] # recorded but NOT frozen — documented drift surface
  artifact_dir, era5_ic_dirs, julia_binary_path
[env_vars] CALIBRATION_CONFIG, CLIMACOMMS_DEVICE, CLIMACOMMS_CONTEXT, ...
```

## ExperimentSpec (config consolidation)

Config files stay executable `.jl` (priors need code), but instead of ~10 loose globals
each must define **one** `const SPEC = CalibrationTools.ExperimentSpec(...)`:

- `name::String`
- `calibrate::CalibrateConfig` (existing struct; its `output_dir` field is now set by
  the launcher, not the config)
- `priors::Vector{ParameterDistribution}`
- `covariance_date_ranges`
- `noise::NoiseSpec` — model_error_scale, regularization quantile, use_latitude_weights,
  min_cosd_lat, rank (currently hard-coded in generate_observations.jl)
- `eki::EKISpec` — scheduler Δt, impose_prior, rng_seed (currently hard-coded in
  run_calibration.jl)
- `preprocessing::PreprocSpec` — pressure_levels, altitude_levels
- `resources::ResourceSpec` — n_workers (derived from priors: 2p+1), walltime, queue,
  device (currently hand-edited in each PBS script)

Backward compat: a `spec_from_legacy_globals()` shim wraps the old-style globals so the
8 existing config files keep working; migrate `lwp_pr.jl` as the template, others lazily.

## New/modified files

**New**
- `src/CalibrationTools.jl` additions (same file, new section) — `ExperimentSpec` +
  sub-structs; `create_snapshot(run_dir; repo_root)` (git archive + dirty patch +
  provenance write); `write_provenance` / `verify_provenance(run_dir)` (rehash
  Manifest/spec/obs files, compare to provenance.toml); `sha256_file` helper.
- `experiments/calibration/amip/launch_calibration.jl` — the single entry point:
  - `julia launch_calibration.jl <config.jl>` → create run dir, snapshot, generate
    observations *as a subprocess with `--project=<snapshot>/experiments/AMIP`*, write
    provenance + launch.pbs, `qsub`.
  - `julia launch_calibration.jl --resume <run_dir>` → `verify_provenance` (hard-fail on
    snapshot tamper), re-`qsub` launch.pbs (ClimaCalibrate resumes from
    last_completed_iteration).
  - `--mode sweep` → same packaging, driver is run_parameter_sweep.jl.
  - launch.pbs template: module load climacommon, campaign julia 1.11.3,
    `JULIA_PROJECT=<snapshot>/experiments/AMIP`,
    `CALIBRATION_CONFIG=<snapshot>/experiments/calibration/amip/config/<cfg>.jl`,
    `CALIBRATION_RUN_DIR=<run_dir>`; deg0025 broken-node guard baked in.

**Modified**
- `experiments/calibration/amip/run_calibration.jl` — read `ENV["CALIBRATION_RUN_DIR"]`
  as output_dir override; take scheduler/rng from `SPEC.eki` instead of hard-coding;
  at startup call `verify_provenance` when running inside a package (presence of
  `package/provenance.toml`) and log the provenance summary into driver.log.
- `experiments/calibration/amip/generate_observations.jl` — take noise params from
  `SPEC.noise`; record obs hashes into provenance.
- `experiments/calibration/amip/config/lwp_pr.jl` — migrate to `SPEC` form (template).
- `experiments/calibration/amip/calibration_report.jl` — accept run_dir alone (spec +
  provenance read from `package/`); add a **Provenance** section (git SHA, dirty flag,
  ClimaAtmos/ClimaCore versions, obs hashes) — this alone would have caught the
  0.41→0.42 drift as a one-line diff between two reports.
- `experiments/calibration/amip/run_parameter_sweep.jl` — spec-driven; runs packaged.
- `experiments/calibration/amip/README.md` — lifecycle doc (launch/monitor/resume/report).

## What still CAN drift (recorded, not frozen — documented in README + provenance)

- ClimaArtifacts campaign dir contents (path recorded; effectively append-only).
- `era5_initial_condition_dir` scratch inputs (paths + file listing hash recorded;
  warning at resume if listing changed).
- The campaign julia binary (path + version recorded).
- Depot package dirs could be removed by `Pkg.gc()` — mitigated by running
  `Pkg.instantiate()` inside the snapshot env at launch (verifies all pinned trees
  present while the originals still exist).

## Phasing (keep diff surface safe; a run may be in flight)

1. **Phase 1 — freeze + provenance (the incident-killer)**: snapshot machinery,
   provenance capture/verify, launcher with legacy-config shim, run_calibration.jl
   output-dir/verify hooks. Old ENV-based workflow keeps working untouched.
2. **Phase 2 — spec consolidation**: ExperimentSpec, migrate lwp_pr.jl, move EKI/noise
   hyperparams out of scripts.
3. **Phase 3 — lifecycle integration**: report provenance section, sweep mode, README.

## Verification

1. **Package a pipeline test**: `launch_calibration.jl config/pipeline_test.jl`
   (TEST_CALIBRATION path already exists; short run). Assert: run dir created;
   `package/code` present and 18 MB-scale; provenance.toml fields populated;
   observation_vec.jld2 generated by the snapshot subprocess (driver.log shows
   snapshot paths).
2. **The rebase-immunity test (the point of it all)**: launch a 2-iteration pipeline
   test; after iteration 1 completes, mutate the working repo (edit a source file that
   changes logged output, or `git checkout` an older commit); `--resume`; assert the
   resumed iteration ran from `package/code` (grep driver.log for the snapshot path)
   and the parameter trajectory equals an uninterrupted control run with the same
   rng_seed — with fixed `rng_seed`, TransformUnscented sigma points are
   prior-determined, so `iteration_001/parameters.toml` should be **byte-identical**.
   Restore the repo afterwards.
3. **Tamper test**: modify a file inside `package/` after launch; `--resume` must
   refuse with a provenance mismatch naming the file.
4. **Worker-env test**: confirm worker logs show `--project=<snapshot>` /
   JULIA_PROJECT=snapshot (workers inherit via qsub -V + JULIA_PROJECT propagation).

## Addendum — design-review findings (second-pass verification)

Sharp details a deeper pass verified; the implementation must respect all of these:

1. **cwd footgun (`src/Input.jl:488`)**: `coupler_toml` entries resolve `isfile(file)`
   against the **cwd first**, falling back to `pkgdir(ClimaCoupler)`. Today's PBS
   scripts `cd` into the live repo, so a snapshot run with a live-repo cwd would
   silently read the LIVE `toml/wxquest_progedmf.toml`. The generated driver.pbs must
   `cd {run_root}/package/ClimaCoupler.jl`.
2. **Working tree is ~38 GB** (untracked outputs/junk); tracked content ~18 MB. So the
   snapshot must be `git archive <sha>` + apply a `git diff HEAD` dirty patch — never
   rsync. **git worktree is also rejected**: it shares `.git` with the live repo, so a
   later gc/branch-delete/rebase can invalidate it, violating freeze semantics, and it
   cannot represent a dirty tree.
3. **Untracked config files are the classic footgun**: `git archive` never captures
   them. `validate_snapshot` must assert every spec-referenced file (config .jl,
   coupler YAML, atmos YAML, each `coupler_toml`, era5 dirs) exists *inside the
   snapshot* at launch — fail in seconds, not hours into a PBS job. Untracked files
   under experiments/config/toml get listed in provenance as `untracked_ignored`.
4. **Output symlink**: config `output_dir`s compute to `joinpath(pkgdir(ClimaCoupler),
   "amip_calibration_<x>")`, which lands *inside the snapshot*. Phase 1: launcher
   creates that path as a symlink → `run_root/output`. Phase 2: spec's explicit
   `run_root` removes the hack.
5. **spec.toml as the frozen artifact** (refines the ExperimentSpec decision): the
   frozen definition should be *data* — diffable, hashable, loadable by the report
   without executing code. `ExperimentSpec(path)` dispatches on extension: `.toml`
   parses directly; `.jl` (all 8 existing configs) is included and its globals
   harvested via a shim, then `write_spec_toml` always emits the canonical
   `package/spec.toml`. Priors serialize as `{name, mean, std, lower, upper}`
   (constrained_gaussian args — exactly what the configs pass).
6. **Freeze assertion, greppable per job**: driver.pbs runs
   `@assert pkgdir(ClimaCoupler) == "<run_root>/package/ClimaCoupler.jl"` before the
   calibration — every driver.log then carries proof the snapshot was used.
7. **Depot-GC guard**: `verify_package` checks the Manifest-pinned package slugs still
   exist under the depot (`ispath(~/.julia/packages/<Pkg>/<slug>)`) and warns —
   `Pkg.gc()` is the one way a frozen env can rot.
8. **Observations are frozen too**: generated once at launch from the snapshot, sha256
   recorded; resume verifies the hash and never regenerates (a resumed run must score
   against the identical observation vector). Each fresh run_root also kills the
   stale-`normalization_stats.jld2` hazard by construction.
9. **driver.pbs exports env explicitly** (JULIA_PROJECT, CALIBRATION_CONFIG,
   CLIMACOMMS_*) rather than trusting `qsub -V` inheritance; provenance snapshots the
   launch ENV for forensics. Side benefit: worker `.julia-*.out` files land in the
   snapshot cwd instead of littering the live repo root.
10. **run_root is immovable**: the serialized interface and driver contain absolute
    run_root paths; moving a package breaks it (document; regeneration on move is out
    of scope).
11. **Resume must handle partially-run members (restart-diagnostics bug, observed
    2026-07-28)**: a member resumed from a mid-run model checkpoint restarts its
    monthly-diagnostics windows at the checkpoint time, so the new `output_active`
    generation holds wrong-dated averages (e.g. a single "monthly" slice at
    2007-11-03 instead of Sep/Oct/Nov-01 means) → observation_map crashes with
    "Too many NaNs", and every naive resubmission re-enters the wedge (the SON
    spin-up run burned 4 driver attempts this way). `--resume` must detect members
    whose checkpoint says "started" with a mid-run model checkpoint and either
    delete their state (rerun from scratch) or verify their diagnostics cover the
    full sample window before proceeding. Sizing driver walltime so iterations
    never straddle job boundaries avoids the case but must not be relied on.
