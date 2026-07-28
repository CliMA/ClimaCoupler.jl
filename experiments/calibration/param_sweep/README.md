# Parameter sweep utility

Run an AMIP parameter sweep from two inputs: a coupler YAML and the
parameters to vary. Parameters mix freely: exact value lists, priors
(swept over their quantiles), and pinned scalars, expanded to the full
factorial. The sweep is a single ClimaCalibrate iteration with no
ensemble update, so each member is an independent forward run of the same
`CouplerModelInterface` that the `amip` calibration uses. The analysis
compares the members' monthly-mean output against near-identical
replicate runs with a latitude-weighted RMS difference, so it measures
whether each parameter moves the model above its internal variability.
There is no comparison to observations.

## Quickstart

1. Clone ClimaCoupler on Derecho.
2. Copy `example_sweep.jl` and edit the inputs. Point `output_dir` at
   scratch.
3. From the repository root:

```
module purge && module load climacommon/2025_02_25
JULIA=/glade/campaign/univ/ucit0011/software/julia/julia-1.11.3/bin/julia
$JULIA --project=experiments/AMIP -e 'using Pkg; Pkg.instantiate()'
$JULIA --project=experiments/AMIP \
    experiments/calibration/param_sweep/run_sweep.jl my_sweep.jl
```

This checks every parameter name against the ClimaParams registry (a
typo is a hard failure before any GPU time), then `add_workers` submits
one GPU worker per member and the `WorkerBackend` runs the members on the
pool as the workers join. Set `n_workers` in the experiment file to use
fewer workers than members; the members then run in waves. The driver has
to stay alive while the members run, so start it in a session that
survives, for example with `nohup` or `screen`. Workers are cancelled when
the driver exits. Pass `--no-submit` to only write the member inputs out.

Rerun the analysis anytime, for example after resubmitting a failed
member. It takes the same experiment file, which is the only record of the
design, so keep it next to the sweep:

```
$JULIA --project=experiments/AMIP \
    experiments/calibration/param_sweep/run_sweep.jl my_sweep.jl --analyze
```

Rerunning the experiment file skips members that ClimaCalibrate recorded
as completed, so a partially failed sweep is resumable with the same
command.

## What a sweep produces

The analysis is a log report. Per variable, it gives the largest RMS
response across the members against member 1 (the base point) and the
replicate noise floor. A signal-to-noise ratio well above 1 means the
parameter moves this variable above internal variability; near or below 1
means its effect is buried. A sweep member that reproduces member 1 bit
for bit is flagged loudly: the parameters it changes are registered but
not wired into the loaded model configuration.

Everything else lands in `output_dir`:

- `iteration_000/member_NNN/parameters.toml`: the values each member ran
  with.
- `iteration_000/member_NNN/`: each member's model output. The forward
  model logs go to the workers instead, in `worker_N.log` next to where
  the driver runs.

## Costs

One member is a 38-day AMIP run (one week of spin-up plus the sampled
month), about 1.2 wall hours on one Derecho GPU node. A sweep runs
`members + replicates` forward models once. Derecho charges exclusive jobs
for all 4 GPUs of a node, so consider adding `q = "preempt"` to
`worker_options` (20% charge): sweep members are independent and
resumable, the workload preemption punishes least.
