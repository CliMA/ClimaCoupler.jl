# AMIP pipeline plan (three-input entry point)

Status: IMPLEMENTED (pipeline.jl, run_amip.jl, example_experiment.toml,
README.md). The implementation deviates from the plan below in one
architectural point: instead of functionizing the three scripts (PR 1),
the pipeline WRITES a self-contained config file into the run directory
and drives the existing scripts through ENV["CALIBRATION_CONFIG"]. That
reuses the battle-tested script paths unchanged and the generated config
doubles as a rerunnable record of the inputs. Deviations from the
interface below: `observations` is required in v1 (every supported
observable has a catalog source, so the response-only sweep is deferred),
and sweep members are given as explicit Dicts (no factorial grids yet).

## Goal

One AMIP-specific entry point that takes three inputs and does the rest:

1. A coupler YAML config (which model to run).
2. Parameters: a prior, or an explicit list of parameter values.
3. Optionally, observations to compare against.

A prior means calibrate. A list means sweep. Observations given as short
names are built automatically from the satellite catalog. Today this
workflow is three scripts, two manual steps, and one edit-this-block
section. The pipeline folds them into one call so that someone outside
this project can run a sweep or a calibration without learning the
internals.

## User interface

```julia
include("experiments/calibration/amip/pipeline.jl")

# Calibration: prior + observations.
amip_pipeline(
    "config/amip_configs/amip_calibration.yml";
    parameters = [
        PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 700, 300, 3600),
        PD.constrained_gaussian("snow_autoconversion_timescale", 1800, 700, 300, 3600),
    ],
    observations = ["lwp", "pr"],
    output_dir = "amip_calibration_myrun",
)

# Sweep: explicit member list, observations optional.
amip_pipeline(
    "config/amip_configs/amip_calibration.yml";
    parameters = [
        Dict("rain_autoconversion_timescale" => 900.0),
        Dict("rain_autoconversion_timescale" => 3600.0),
        Dict("rain_autoconversion_timescale" => 14400.0),
    ],
    replicates = 4,
    observations = ["lwp", "pr"],
    output_dir = "amip_sweep_myrun",
)
```

Dispatch on `parameters`:

- `Vector{ParameterDistribution}` (or one combined distribution):
  calibration. TransformUnscented with 2p + 1 members, the validated EKP
  settings (DefaultScheduler(0.1), impose_prior), packed GPU workers.
  Requires `observations`.
- `Vector{<:AbstractDict}`: sweep. One member per entry, plus `replicates`
  near-default members that measure the internal-variability noise floor.
  With `observations` the report gates signal against noise per variable
  (the identifiability-map machinery). Without `observations` pass
  `short_names` instead, and the report shows response spread only.

`observations` accepts:

- `Vector{String}`: short names resolved through the existing loader
  catalog (lwp, cl, pr, swcre, lwcre, ERA5 fields). The pipeline builds
  `observation_vec.jld2` and the coverage masks automatically, or reuses
  them when already present in `output_dir`.
- A path to a prebuilt `observation_vec.jld2`.
- `nothing` (sweep only).

Keyword defaults carry the validated setup and every one is overridable:
`sample_date_ranges` (Oct 2010), `covariance_date_ranges` (Oct 2006 to
2010), `spinup` (7 days), `extend` (1 month), `n_iterations` (8),
`noise_groups` (0.2 floor, radiation and cl variables 0.5),
`reduction` (`:zonal` or an integer coarsening factor), `rng_seed`,
`workers_per_node` (4), `walltime`.

Every run passes preflight before any worker is submitted. A preflight
failure aborts with the check table. `skip_preflight = true` exists for
debugging.

## What exists and what is new

All five scripts already separate reusable functions from an
`if abspath(PROGRAM_FILE)` entry block, so the work is moving script-body
logic into callable functions, then one dispatch layer.

**PR 1, mechanical functionization (no behavior change):**

- `generate_observations.jl`: script body becomes
  `build_observation_vector(cfg; noise_groups, covariance_date_ranges)`.
  The script entry calls it.
- `run_calibration.jl`: EKP construction, worker startup, and `calibrate`
  become `calibrate_amip(cfg, priors; eki_kwargs...)`.
- `run_parameter_sweep.jl`: the edit-this-block sweep specification
  becomes the argument of `sweep_amip(cfg, points; replicates,
  observations)`. Commit the script first, it is currently untracked.
- Existing config files and PBS scripts keep working unchanged.

**PR 2, the pipeline:**

- `pipeline.jl`: `amip_pipeline(yaml; parameters, observations, ...)`,
  which builds the `CalibrateConfig`, materializes observations, runs
  preflight, dispatches on the parameter type, and writes a run summary
  (parameters, observations, settings, git SHA) into `output_dir`.
- `run_amip.jl`: thin driver so the pipeline also works from a PBS job:
  `julia run_amip.jl experiment.toml`, where the TOML carries the three
  inputs for people who prefer a file to a REPL call. The TOML supports
  priors as `[name, mean, sigma, lower, upper]` rows.
- A PBS template generator for the CPU-driver launch pattern, so `qsub`
  is the only scheduler command a user touches.
- README quickstart with the two examples above.

## Relation to existing plans

This is the slim version of the packaging plan's Phase 2
(ExperimentSpec): it consolidates the knobs into one call instead of one
struct, and does not freeze code. The on-hold packaging framework can
later wrap `run_amip.jl` as its driver without rework. The sweep half is
also the place where the generic core could move upstream into
ClimaCalibrate later; the pipeline keeps that split clean by calling
`sweep_amip` rather than inlining it.

## Verification

1. Pipeline test config (TEST_CALIBRATION path) through both modes:
   a 2-point sweep and a 1-parameter calibration on emulated diagnostics.
2. Reproduce the lwp+pr calibration through `amip_pipeline` with
   `rng_seed = 42` and assert the iteration-1 parameter files match the
   existing run.
3. A sweep with `observations = nothing` and `short_names = ["lwp"]`
   produces the response-only report.
4. Preflight failure (a not-wired parameter) aborts before any qsub.
