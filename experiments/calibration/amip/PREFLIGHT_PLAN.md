# Pre-flight validator plan (v1)

Status: IMPLEMENTED. Validated on Derecho (job 6936450): all checks pass on config/phase2.jl, and the wiring check fails on the known not-wired parameter mixing_length_diss_coeff.

## Goal

Catch configuration mistakes before a calibration spends a 12 hour GPU
allocation. Every check below corresponds to a failure we hit by hand:

- `mixing_length_diss_coeff` produced members bit-identical to the center
  because the parameter is not wired in Atmos 0.42. We lost a full sweep.
- The `cl` coverage mask silently applied to zero points because the
  observation dimension was `height` and the simulation dimension was `z`.
- Observation vectors were sanity-checked by hand after each generation
  (length, per-variable noise level, mask sizes).
- Worker count, walltime, and ensemble size were checked by hand against
  the priors and the measured model speed.
- Resuming into a directory with partial members triggers the restart
  diagnostics bug (wrong-dated monthly means, then a NaN crash loop).

## Deliverable

One new script and no changes to the calibration driver:

```
experiments/calibration/amip/preflight.jl
```

Usage, from the repo root:

```
julia --project=experiments/AMIP experiments/calibration/amip/preflight.jl \
    experiments/calibration/amip/config/phase2.jl
```

The script includes the config file, runs all checks, and prints one line
per check: `PASS`, `WARN`, or `FAIL` with a short reason. It exits with a
nonzero code if any check fails. The launch procedure becomes: generate
observations, run preflight, submit only on success.

## Checks

### A. Config consistency (static, runs in seconds)

1. Priors: each prior mean is inside its bounds, names are unique, and
   names are nonempty.
2. Ensemble and resources: for `TransformUnscented` the ensemble size is
   `2p + 1`. Print the required worker count. Estimate the wall time per
   iteration from the simulation length and the measured speed (about
   2.0 simulated years per day, overridable with a flag). `WARN` if the
   estimate for all iterations exceeds the 720 minute worker walltime in
   `run_calibration.jl`.
3. Dates: every `sample_date_range` and covariance date range lies inside
   the observation record of each data loader. At least 2 covariance
   samples exist.
4. Grid: if the config sets `COARSEN_FACTOR`, confirm it divides the
   comparison grid in both directions.

### B. Observation vector (reads the generated files, runs in seconds)

Runs only if `observation_vec.jld2` exists, otherwise prints `WARN` with
the generation command.

1. The number of samples equals `length(sample_date_ranges)`.
2. The constraint blocks match `short_names` (names and count), and the
   total length matches the expected grid size per variable.
3. No block is all NaN and the covariance diagonal is finite and positive.
4. Per-variable relative noise (sigma over signal RMS) is printed. `WARN`
   outside 5 to 200 percent.
5. Each coverage mask masks at least one point and not all points. The
   mask dimensions match the simulation dimensions after canonical name
   matching, so a `height` versus `z` mismatch is a `FAIL`, not a silent
   skip.

### C. Run directory hygiene (runs in seconds)

1. If `output_dir` already contains iterations, find partial members
   (member directories without complete diagnostics). Partial members are
   a `FAIL` with the cleanup instruction, because resuming over them
   triggers the restart diagnostics bug.
2. If `ENV["CALIBRATION_CONFIG"]` is set, it must point at the config
   under check.

### D. Parameter wiring (constructs parameter structs, runs in minutes)

For each prior, confirm the parameter actually reaches the model:

1. Write two member TOML files: the prior center, and the center with
   this one parameter moved by one constrained sigma (clipped inside the
   bounds).
2. Build the merged parameter dictionary exactly the way a member does:
   `Input.get_coupler_config_dict`, then
   `CalibrationTools.add_parameter_filepath!`, then the same
   `merge_toml_files` and `merge_override_default_values` path used in
   `ext/ClimaCouplerClimaAtmosExt.jl`.
3. Construct the ClimaAtmos parameter struct from each merged dictionary
   with the model choices from the config. Do not build spaces and do not
   call `get_simulation`, so no GPU is needed.
4. Recursively diff the two structs over all numeric fields. No
   difference means the parameter never lands in the model, which is a
   `FAIL`.

Validation requirement for this check: it must `FAIL` on
`mixing_length_diss_coeff` under Atmos 0.42 (the known not-wired case,
where the dissipation coefficient is derived from other parameters) and
`PASS` on every parameter in `config/phase2.jl`. If the struct-level diff
does not catch the known case, fall back to a short-run diff (a few model
steps, hash of the prognostic state) as the check, and note the cost.

Escape hatch: a parameter list `PREFLIGHT_NON_ATMOS_PARAMS` in the config
for land or coupler parameters that the atmos struct check cannot see.
These print `WARN` (unchecked), not `FAIL`.

## Implementation notes

- The script reuses the config include mechanism from
  `run_calibration.jl`, so a config needs no changes to be checkable.
- Checks A to C need only ClimaAnalysis, JLD2, and the config. Check D
  loads ClimaAtmos through `code_loading.jl` and is the slow part. A
  `--skip-wiring` flag runs A to C alone for quick re-checks.
- The recursive struct diff and the partial-member finder go into
  `src/CalibrationTools.jl` so tests can call them. Everything else stays
  in `preflight.jl`.
- Output is log only, matching the steering indicators decision. No files
  are written except the two temporary TOML files in a scratch directory.

## Tests

1. Unit tests for the struct diff and the partial-member finder.
2. `preflight.jl config/pipeline_test.jl` passes end to end.
3. Negative controls: a config with a prior named
   `mixing_length_diss_coeff` fails check D. A truncated observation
   vector fails check B. A run directory with a fabricated partial member
   fails check C.

## Out of scope (later PRs)

- The reusable run watcher and failure-signature monitor.
- `steering_history`: extract and tabulate steering blocks from a
  driver log.
- `sweep_report`: standard admit or reject gates over a sweep directory.
- The packaging framework (on hold) absorbs launch, provenance, and
  resume. Preflight stays independent so it can run with or without
  packaging.
