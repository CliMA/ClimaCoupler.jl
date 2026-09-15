# Parameter sweep utility

Run an AMIP parameter sweep from two inputs: a coupler YAML and the
parameters to vary. Parameters mix freely: exact value lists, priors
(swept over their quantiles), and pinned scalars. The members are every
combination of the swept values. The sweep is a single ClimaCalibrate
iteration with no ensemble update, so each member is an independent
forward run of the same `CouplerModelInterface` that the `amip`
calibration uses. Each member writes its own model output; comparing the
members is left to you.

A prior is swept over evenly spaced quantiles, `i / (prior_points + 1)`
for `i in 1:prior_points`, so the default of 3 points gives the 25th,
50th, and 75th percentiles. The endpoints are left out because a prior's
0th and 100th percentiles are often unphysical.

## Quickstart

1. Clone ClimaCoupler on Derecho.
2. Copy `example_sweep.jl` and edit the inputs. Point `output_dir` at
   scratch.
3. From the repository root:

```
module purge && module load climacommon/2025_02_25
julia --project=experiments/AMIP -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/AMIP \
    experiments/calibration/param_sweep/run_sweep.jl my_sweep.jl
```

This checks every parameter name against the ClimaParams registry (a
typo is a hard failure before any GPU time), then `add_workers` submits
one GPU worker per member and the `WorkerBackend` hands members to
workers as they free up. The driver has to stay alive while the members
run, so start it in a session that survives, for example with `nohup`,
`screen`, or `tmux`. Workers are cancelled when the driver exits. Pass
`--no-submit` to only write the member inputs out.

Set `n_workers` to run the members on fewer workers than there are
members. Pick a divisor of the member count: any remainder leaves a final
wave with idle workers, and a wave takes as long as its slowest member
either way. With 9 members, 9 or 3 workers keep every worker busy, while
8 leaves one member to run alone in a second wave and roughly doubles the
wall time.

Rerunning the experiment file skips members that ClimaCalibrate recorded
as completed, and a member that stopped partway resumes from its last
checkpoint, so a partially failed sweep is resumable with the same
command. A resumed member's diagnostics are affected by the restart, see
Costs below. The experiment file is the only record of the design, so
keep it next to the sweep.

## What a sweep produces

Everything lands in `output_dir`:

- `iteration_000/member_NNN/parameters.toml`: the values each member ran
  with.
- `iteration_000/member_NNN/`: each member's model output, holding the
  diagnostics the coupler YAML asks for. The forward model logs go to the
  workers instead, in `worker_N.log` next to where the driver runs.

## Costs

One member is a 38-day AMIP run (one week of spin-up plus the sampled
month), about 1.2 wall hours on one Derecho GPU node. The member count is
capped at 64 to catch a value list that accidentally asks for hundreds of
these; raise `max_members` if a large sweep is what you want.

Derecho charges exclusive jobs for all 4 GPUs of a node, so `q = "preempt"`
in `worker_options` costs 20% of the regular charge, at the cost of an
unpredictable walltime.

A preempted member loses at most 10 simulated days of stepping, but the
diagnostics are not that safe. Checkpoints are every 10 simulated days
and do not line up with the monthly averaging windows the diagnostics
use, and a checkpoint does not store the partially accumulated averages.
A month that a restart falls inside is therefore averaged over only the
part after the restart. Use `preempt` when the sweep is exploratory, and
the regular queue when the diagnostics have to be comparable across
members.
