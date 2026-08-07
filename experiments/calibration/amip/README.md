# AMIP calibration

Calibrate the coupled AMIP model against satellite observations on
Derecho.

## Start here

1. This README: how to launch a run and operate it on Derecho.
2. `CALIBRATION_GUIDE.md`: the methodology. Gates before a run, the
   calibration ladder, how to read a running calibration.

The campaign's run-by-run configs are in `config/`, each header holding
that run's design, predictions, and verdict. A retired config is
evidence, not dead code.

## Two ways to launch

### 1. The pipeline (new experiments)

One TOML with three inputs: a coupler YAML, parameters, and observations.

1. Copy `example_experiment.toml` and edit the three inputs. Point
   `output_dir` at scratch.
2. From the repository root, on a `develop`-queue node:

```
module purge && module load climacommon/2025_02_25
JULIA=/glade/campaign/univ/ucit0011/software/julia/julia-1.11.3/bin/julia
$JULIA --project=experiments/AMIP -e 'using Pkg; Pkg.instantiate()'
$JULIA --project=experiments/AMIP \
    experiments/calibration/amip/run_amip.jl my_experiment.toml
```

The pipeline writes a self-contained config into `output_dir`, builds the
observation vector, runs a preflight check, and submits the PBS driver.
Pass `--no-submit` to prepare everything and print the `qsub` command.

One `[[priors]]` table per free parameter (TransformUnscented,
`2p + 1` members for `p` parameters). From Julia, the same interface is
`amip_pipeline` in `pipeline.jl`.

### 2. The campaign path (hand-written config)

Every campaign run uses a config file in `config/` selected through
`ENV["CALIBRATION_CONFIG"]`, in four steps:

1. Write the config: `config/<name>.jl`, a coupler YAML in
   `config/amip_configs/`, and a fixed-parameter TOML in `toml/` if
   parameters move between fixed and free. Use the most recent config in
   `config/` as the template.
2. Create the run dir on scratch and symlink it from the repository
   root. This is mandatory: a real run dir in the repository fills the
   100 GiB home quota and truncates JLD2 checkpoints.
3. Prep job (develop queue, 1 GPU): run `generate_observations.jl`, then
   `preflight.jl <config>`. Submit only on green. The wiring check
   matters most: a parameter name the model silently ignores must be a
   hard failure now, not a flat posterior six iterations later.
4. Start the driver: `run_calibration.jl` under
   `--project=experiments/AMIP` with `CALIBRATION_CONFIG` set. Then,
   before the members burn GPU time, read
   `iteration_001/member_001/parameters.toml` and the priors plot and
   confirm they match the config.

After iteration 1 completes, check the go/no-go gate: each observable's
iteration-1 ensemble spread of G must beat its noise floor, or the
calibration cannot move that observable.

## What a run produces

Everything lands in `output_dir`:

- `observation_vec.jld2`, `coverage_masks.jld2`: the observation.
- `driver.log`: driver progress, with per-iteration steering indicators.
- `iteration_NNN/`: member output, `eki_file.jld2`, `G_ensemble.jld2`,
  per-iteration `g_vs_obs.png` and bias plots, evolution GIFs.

Reports: run `calibration_report.jl` as a PBS analysis job (develop
queue, `select=1:ncpus=4:ngpus=1:mem=40GB`), not on a login node.

## Derecho operations

- Queues: `main` charges exclusive nodes for all 4 GPUs regardless of
  the request; `develop` (6 h limit) is shared. Workers pack 4 per GPU
  node, which cuts charged GPU-hours 4x.
- Cost scale: one member is a 38-day AMIP run, about 1.2 GPU hours. A
  calibration runs `(2p + 1) * n_iterations` members.
- The driver is cheap and can run in tmux on a login node instead of a
  CPU job. Three constraints:
  - Workers must bind to the HSN network. Set
    `CALIBRATION_WORKER_EXENAME=julia_worker_hsn_bind.sh` (in this
    directory). Without it every worker fails with "Master process
    could not connect within 60.0 seconds".
  - Login nodes enforce a 10 GiB per-user memory cgroup. Start the
    driver with `--heap-size-hint=3G` or the OOM killer takes it.
  - Julia buffers log output for hours. Monitor by process liveness and
    member file-write recency, never by log tail alone.
- A calibration is resumable: rerun the driver with the same config and
  it continues from the last completed iteration. Clean partial members
  first; a member resumed from a mid-run model checkpoint writes
  wrong-dated monthly diagnostics and crashes the observation map.
- Login-node `qstat` is a cache. It can report "Unknown Job Id" for a
  job you just submitted and freeze counters for 30 minutes. Wait on
  the job's own log markers instead.
