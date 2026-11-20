# Precompilation workloads

Building a `CoupledSimulation` and taking the first `step!` takes tens of minutes, almost all
of it single-core compilation repeated on every launch. These packages run a short,
production-configured coupled simulation while they precompile, so that code is saved in
package images and a run starts in minutes. Building the images takes ~95 min, once.

**Use it** for repeated launches of unedited code: CI, long runs, production, calibration.
**Don't** while editing ClimaCoupler, ClimaAtmos or ClimaCore: every edit invalidates the
images and triggers a ~95 min rebuild. That is why this is a separate environment, not a
dependency of `experiments/CMIP` or `experiments/AMIP`.

## Usage

```bash
# 1. Build. CUDA is required (otherwise the images hold CPU code); one task at a time,
#    since each build holds a ~9 GB simulation on the GPU.
CLIMACOMMS_DEVICE=CUDA JULIA_NUM_PRECOMPILE_TASKS=1 \
  julia --project=experiments/warmup -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# 2. Check. Always run: a failed workload does not fail Pkg.precompile().
julia --project=experiments/warmup experiments/warmup/check_warmup.jl

# 3. Run, loading CMIPWarmup (and with it every image) first.
CLIMACOMMS_DEVICE=CUDA julia --project=experiments/warmup \
  -e 'using CMIPWarmup; include("experiments/CMIP/run_simulation.jl")' \
  -- --config_file <your config>.yml
```

On a GPU too small for one coupled simulation, add `WARMUP_H_ELEM=8` to step 1. `h_elem`
only sets array sizes, so these images still serve production runs at `h_elem = 16`.

## What invalidates the images

A run that compiles different code silently recompiles from scratch.

* **Model configuration.** `configs/*.yml` mirror production; model choices (float type,
  topography, gravity-wave model, ...) are part of the compiled types.
* **Vertical resolution.** `z_elem` is a type parameter (`Nv` in `VIJFH{FT, Nv, ...}`), so
  there is no `z_elem` override.
* **Package versions.** `Manifest-v1.13.toml` is `experiments/CMIP/Manifest-v1.13.toml` plus
  the warmup packages. When the CMIP manifest changes, re-seed and rebuild:

  ```bash
  cp experiments/CMIP/Manifest-v1.13.toml experiments/warmup/Manifest-v1.13.toml
  julia --project=experiments/warmup -e 'using Pkg; Pkg.resolve()'
  ```

`check_config_drift.jl` and `check_version_drift.jl` catch the first and last.

## Packages

| package | caches | depends on |
|---|---|---|
| `CouplerWarmupCommon` | shared setup for the component workloads | ClimaCoupler |
| `AtmosWarmupSetup` | `ClimaAtmos.get_simulation` | — |
| `AtmosWarmupExplicit` | explicit tendencies (`remaining_tendency!`) | Setup |
| `AtmosWarmupStep` | implicit solve, Jacobian, timestepper | Setup, Explicit |
| `OceanWarmup` | the Oceananigans ocean as the coupler builds it | Setup |
| `LandWarmup` | ClimaLand as the coupler builds it | Setup |
| `SeaIceWarmup` | ClimaSeaIce (needs a built ocean) | Setup, Ocean |
| `CMIPWarmup` | the coupled `run!`: field exchange, fluxes, exchange grid | all of the above |

A package image cannot exceed ~2 GiB, hence the split. Edges follow the call graph so each
package reuses the code cached beneath it; removing one pushes a package over the limit.
Every output and physics schedule is shortened to fire at least once, since code behind a
schedule that never fires is never compiled.

## Troubleshooting

`check_warmup.jl` reports a workload that threw (recorded in `workload_report.txt`), an image
over the size limit (stale on disk despite `status = ok`), and an image that cached almost
nothing. Fix the lowest failing package in the table first, since the packages above it then
cache its code themselves. Delete the affected images before rebuilding: Pkg considers an
image from a failed workload up to date.

Timings, image sizes and limits are in [`precompile_info.md`](precompile_info.md).
