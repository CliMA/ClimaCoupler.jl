# Calibration run registry plan (ClimaCalibrate)

Status: PLANNED. Target: one PR against ClimaCalibrate main (a4072867, which
contains the merged async-worker machinery we run on).

## Goal

A database of calibration runs, so that finding and tracking runs needs no
hand-written monitors. Today the mapping between a run directory, its PBS
jobs, its config, and its progress lives only in my head and in qstat. The
registry stores that mapping and derives the live status on demand.

The user-facing result is one call:

```
julia> ClimaCalibrate.Registry.status()
name                     state    iter   driver          workers  updated
amip_calibration_phase2  RUNNING  3/8    6935308 (R)     9 alive  4 min ago
amip_calibration_lwp_pr  DONE     10/10  finished        -        21 days ago
amip_calibration_lwp_cl  DEAD     2/10   6912011 (gone)  0        30 days ago
```

`DEAD` means the ledger has no final event and the driver is not alive. That
is the state that today requires detective work to notice.

## Design

Two plain-text pieces, both TOML, no new dependencies (TOML and Dates are
already dependencies):

1. A per-run ledger `<output_dir>/calibration_ledger.toml`. Append-only
   `[[event]]` blocks. Appending a block keeps the file valid TOML, and
   `TOML.parsefile` returns the event list. The ledger travels with the run.
2. A central index, default `~/.climacalibrate/runs.index`, one absolute
   run directory per line. Written once per registration with an atomic
   append. It exists only so `status()` can find runs without scanning
   scratch. The home directory is shared across Derecho nodes, so every
   login and compute node sees the same index. `ENV["CLIMACALIBRATE_REGISTRY"]`
   overrides the location (for a lab-shared index).
   `ENV["CLIMACALIBRATE_NO_REGISTRY"] = "1"` disables all writes.

### Events

Every event has `event`, `time`, `host`, `pid`. Specific events add:

- `registered`: output_dir, n_iterations, ensemble_size, backend type,
  driver scheduler job id (from `PBS_JOBID` or `SLURM_JOB_ID`), Julia
  project path, and the worker job name `julia-<pid>`. A resumed run appends
  a second `registered` event to the same ledger. The index deduplicates by
  `realpath`.
- `workers_added`: count, device, and the scheduler job ids (PBS via the
  existing `_pbs_worker_job_ids(worker_jobname())`, Slurm via
  `squeue --name`).
- `iteration_completed`: iteration number. Written after
  `observation_map_and_update!` succeeds, so it means the EKP update is on
  disk, not just the forward models.
- `finished`: normal loop exit.
- `terminated`: the EKP termination condition fired early.
- `failed`: exception text, appended from a catch block before rethrow.

### Status derivation (read side, no daemon)

`status()` reads the index, then each ledger, newest event last:

- final event present: `DONE`, `TERMINATED`, or `FAILED`.
- no final event: query the driver job id with `qstat`/`squeue`. Alive means
  `RUNNING` (plus `Q` shows as `QUEUED`). Gone means `DEAD`. When the run
  has no scheduler job id (interactive driver), fall back to
  `kill -0 pid` on the recorded host name matching the current host,
  otherwise report `UNKNOWN`.
- iteration progress comes from the ledger, cross-checked against
  `last_completed_iteration(output_dir)` on disk (the disk wins, since it is
  what resume uses).
- live worker count via `qselect -N julia-<pid>` with the recorded pid.

Scheduler queries are isolated in one function
`scheduler_job_state(job_id) -> :running | :queued | :gone | :unknown` so
tests can inject a fake.

## Code changes (all in ClimaCalibrate)

New file `src/registry.jl`, module `Registry` (~200 lines):

- `register_run!(output_dir; n_iterations, ensemble_size, backend)`
- `record_event!(output_dir, event; pairs...)` with an flock-free atomic
  append (open with append flag, write one full block, close)
- `annotate!(output_dir, pairs...)` for caller metadata. Our ClimaCoupler
  driver will use it to record the CALIBRATION_CONFIG path and the parameter
  names, so `status()` output is self-describing.
- `runs()`, `status(io = stdout)`, `run_state(output_dir)`
- included from `ClimaCalibrate.jl`, exported as `Registry`

Hooks, each one or two lines:

- `src/calibration.jl`: the three `calibrate` methods share the loop shape,
  so add `register_run!` before the loop, `record_event!(:iteration_completed)`
  after each update, `:finished` or `:terminated` after the loop, and a
  `try ... catch` that records `:failed` and rethrows. To avoid triplicating
  this, extract the shared loop into one internal
  `run_calibration_loop(run_iter, ekp, prior, output_dir, n_iterations, interface)`
  where `run_iter` is a closure over each backend's `run_iteration`
  arguments. The loop bodies are already identical, only the
  `run_iteration` signatures differ (the HPC variant takes the interface
  file path and exe flags).
- `src/backends/workers.jl`: at the end of `_add_workers`, record
  `workers_added`. This runs on the spawned task, which is safe because
  appends are atomic single writes.
- Registry failures must never kill a calibration: every write is wrapped so
  an error becomes one `@warn`.

Not in this PR (later, if wanted): member-level `JobInfo` recording for the
HPC backends, stall detection (running but no iteration progress for N
hours), and a standalone lister that avoids loading the full package.

## Testing

- Ledger round-trip: append events, parse, assert order and fields.
- State machine: fabricate ledgers for each terminal and non-terminal case
  with an injected fake `scheduler_job_state`, assert the derived state.
- Resume: two `registered` events, index stays deduplicated.
- Concurrency: many tasks appending in parallel, parse succeeds, no torn
  blocks.
- Integration: the existing `JuliaBackend` end-to-end test asserts the
  ledger contains `registered`, `iteration_completed` per iteration, and
  `finished`. The worker-backend test asserts `workers_added`.
- Opt-out: with `CLIMACALIBRATE_NO_REGISTRY=1` no files appear.

## Rollout for our runs

1. Branch `ne/run-registry` off ClimaCalibrate main, implement, PR upstream.
2. Point ClimaCoupler's `experiments/AMIP` at the branch (same `[sources]`
   mechanism as ne/async used) once ClimaCoupler moves to a ClimaCalibrate
   release that includes the merged async code.
3. Add `Registry.annotate!` and the config path to our `run_calibration.jl`.
4. Retire the hand-written background monitors: the watcher reduces to
   `Registry.status()` plus the steering log.

Note: our current runs pin ClimaCalibrate to the ne/async branch state, so
the registry only applies to runs launched after the dependency moves. Do
not move the pin while phase2 and the 2D validations are in flight.
