# Precompilation workloads

A coupled run spends most of its startup compiling, not computing: building the
`CoupledSimulation` and taking the first `step!` take tens of minutes, almost all of it
compilation, while every later step takes a fraction of a second. That compilation runs on a
single core and is repeated on every launch.

The packages in this directory run a short, production-configured coupled simulation while
they precompile. The compiled setup and first-step code is saved into package images on
disk, so a run that loads them starts in minutes instead of tens of minutes.

You pay for this once, when building the images (~95 min). Every launch after that is faster,
as long as the code, the package versions and the model configuration stay the same.

## When it is worth it

**Use it** when you launch the same, unedited code many times: CI, long runs, production
runs, and calibration (each persistent worker pays startup only once).

**Do not use it** while you are editing ClimaCoupler, ClimaAtmos or ClimaCore. These packages
are built on top of them, so every edit throws the images away and the next
`Pkg.precompile()` rebuilds all of them (~95 min). This is why the warmup packages live in
their own environment and are *not* dependencies of `experiments/CMIP` or
`experiments/AMIP`: nobody pays the rebuild cost unless they choose to.

## Usage

### 1. Build the images

```bash
CLIMACOMMS_DEVICE=CUDA JULIA_NUM_PRECOMPILE_TASKS=1 \
  julia --project=experiments/warmup -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

* `CLIMACOMMS_DEVICE=CUDA` is required. Without it the images contain CPU code, and a GPU
  run gains nothing.
* `JULIA_NUM_PRECOMPILE_TASKS=1` builds one package at a time. Each build holds a full
  production-resolution simulation (~9 GB) on the GPU, so two at once exceed a 15 GB card.
  Raise it only if your GPU can hold that many simulations.

On a GPU too small for even one coupled simulation, add `WARMUP_H_ELEM=8`; see
[Building on a small GPU](#building-on-a-small-gpu).

### 2. Check the build

```bash
julia --project=experiments/warmup experiments/warmup/check_warmup.jl
```

Always run this. A failed workload does not fail `Pkg.precompile()`, so a build can look
fine and still give you no speedup. The script exits non-zero if any image is unusable; see
[Troubleshooting](#troubleshooting).

### 3. Run with the images loaded

Run in the `experiments/warmup` environment and load `CMIPWarmup` before anything else.
Loading it loads all the images beneath it:

```bash
CLIMACOMMS_DEVICE=CUDA julia --project=experiments/warmup \
  -e 'using CMIPWarmup; include("experiments/CMIP/run_simulation.jl")' \
  -- --config_file <your config>.yml
```

The `experiments/warmup` environment contains everything `experiments/CMIP` does, pinned to
the same versions. Loading the images adds to package-load time, but far less than the
compilation they save.

## What the images can and cannot cover

An image only helps a run that compiles *exactly the same code*. If anything below differs,
the run quietly recompiles from scratch, which looks just like "precompilation didn't help".

**Model configuration.** `configs/atmos_warmup.yml` and `configs/cmip_warmup.yml` mirror the
production configuration. Many model choices (float type, topography, gravity-wave model,
and so on) are part of the compiled types, so a run with different choices cannot use the
images. The images serve the production configuration only. To support another one, the
configs here have to describe it.

**Vertical resolution.** The number of vertical levels (`z_elem`) is a type parameter of
every field (`Nv` in `ClimaCore.DataLayouts.VIJFH{FT, Nv, ...}`). Images built at one
`z_elem` are useless at any other, even one level different.

**Horizontal resolution is free.** `h_elem` only changes array sizes, so images built at one
`h_elem` serve runs at any other.

**Package versions.** The images are tied to the exact package versions they were built
with. `Manifest-v1.13.toml` here is therefore a copy of `experiments/CMIP/Manifest-v1.13.toml`
plus the warmup-only packages. After the CMIP manifest changes, re-seed it and rebuild:

```bash
cp experiments/CMIP/Manifest-v1.13.toml experiments/warmup/Manifest-v1.13.toml
julia --project=experiments/warmup -e 'using Pkg; Pkg.resolve()'
```

Two scripts catch drift before it costs you a slow run:

```bash
julia --project=experiments/warmup experiments/warmup/check_config_drift.jl   # configs vs production
julia --project=experiments/warmup experiments/warmup/check_version_drift.jl  # versions vs experiments/CMIP
```

### Building on a small GPU

The coupled build holds the atmosphere, ocean and sea ice on the GPU at once, which at
production resolution exceeds a 15 GB card. Because horizontal resolution does not affect
the images, you can build at a lower one:

```bash
CLIMACOMMS_DEVICE=CUDA JULIA_NUM_PRECOMPILE_TASKS=1 WARMUP_H_ELEM=8 \
  julia --project=experiments/warmup -e 'using Pkg; Pkg.precompile()'
```

Images built at `h_elem = 8` serve a production run at `h_elem = 16`. There is deliberately
no equivalent for `z_elem`, since lowering it would produce images no production run can
use. The value used is recorded in `workload_report.txt`.

## Package layout

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

The work is split across several packages because a single package image cannot exceed
~2 GiB. The dependency edges follow the call graph (for example, `step!` calls
`remaining_tendency!`, so `AtmosWarmupStep` depends on `AtmosWarmupExplicit`). Each package
therefore reuses the code cached beneath it instead of caching it again. Removing an edge
pushes a package over the limit. As a result the atmosphere packages build one after
another, and only the independent components (ocean, land) can build in parallel.

Each workload runs briefly but at production resolution, with every output and physics
schedule (diagnostics, radiation, cloud fraction, gravity-wave drag, checkpointing)
shortened so it fires at least once. Code behind a schedule that never fires is never
compiled.

## Troubleshooting

`check_warmup.jl` reports three failures that are otherwise invisible:

* **A workload threw.** Each workload catches its own error so that one failure does not
  stop the rest of the build, and records it in `workload_report.txt`.
* **A package image exceeded the ~2 GiB limit.** Linking happens after the workload
  succeeds, so the report says `status = ok` while the image on disk is stale.
* **An image cached almost nothing.** It built without error but is far below its expected
  size.

A failure low in the table makes every package above it cache that code itself, which can
push them over the size limit too. Fix the lowest failing package first. Then delete the
affected images, because Pkg considers an image from a failed workload up to date and will
not rebuild it on its own, and build again.

## Measurements

The timings, image sizes and limits behind these rules are collected in
[`precompile_info.md`](precompile_info.md).
