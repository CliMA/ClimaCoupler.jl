# Precompilation workloads

Packages whose `@compile_workload` runs real model code, so that the inference and codegen
for setup and the first step land in a package image instead of being redone on every launch.

**These are deliberately *not* dependencies of `experiments/CMIP` or `experiments/AMIP`.**
See "Do not put these in a developer's environment" below.

## Why

Nothing in the CliMA stack ships `@compile_workload`s, so each launch re-infers and re-emits
the whole setup and step path on one core. Measured on a coupled CMIP configuration
(prognostic EDMF, Float64, CUDA; see `experiments/CMIP/startup/`):

| phase | wall | of which compilation |
|---|---:|---:|
| package loading | 17 s | — |
| `CoupledSimulation(...)` | 61 min | 84 % |
| first `step!` | 19–23 min | 91 % |
| second `step!` | 0.2 s | 0 % |

Startup is ~99 % single-core: Julia serializes codegen process-wide, so the only parallelism
available is across *processes*, which is what `Pkg.precompile()` provides.

## Layout

| package | caches | depends on |
|---|---|---|
| `CouplerWarmupCommon` | shared prologue for the component workloads | ClimaCoupler |
| `AtmosWarmupSetup` | `ClimaAtmos.get_simulation` | — |
| `AtmosWarmupExplicit` | `remaining_tendency!` | Setup |
| `AtmosWarmupStep` | implicit solve, Jacobian, timestepper | Setup, Explicit |
| `OceanWarmup` | Oceananigans ocean as the coupler builds it | Setup |
| `LandWarmup` | ClimaLand as the coupler builds it | Setup |
| `SeaIceWarmup` | ClimaSeaIce (needs a built ocean) | Setup, Ocean |
| `CMIPWarmup` | the coupled `run!`: field exchange, fluxes, exchange grid | all of the above |

Measured image sizes (production resolution): 1260 / 1205 / 925 / 1240 / 1332 / 108 /
1554 MiB. `SeaIceWarmup` is small *because* `OceanWarmup` is beneath it -- when that build
failed, it grew to 713 MiB re-caching the ocean itself.

## The dependency edges are required, not an optimisation

A package image cannot exceed ~2 GiB. Above that, code in the image references the
serialized blob with 32-bit ELF relocations (`R_X86_64_PC32`, `GOTPC32`) that cannot encode
the displacement. Over the limit `lld` errors out; just under it, the link succeeds silently
and the image segfaults on load inside `jl_restore_package_image_from_stream`.

Measured, same three workloads, three dependency shapes:

| shape | result |
|---|---|
| chain: Setup -> Explicit -> Step | **works**: `.ldata` 1191 / 1193 / 900 MiB, loads |
| fan-out: Setup <- {Explicit, Step} | Step fails: 20 relocation errors |
| fully independent | Step fails: 20 relocation errors |
| single unsplit workload | fails: one 3.05 GiB `sysimg.o`, `lld` refuses |

The edges must mirror the **call graph**. `step!` calls `remaining_tendency!`, so
`AtmosWarmupStep` must depend on `AtmosWarmupExplicit` or it re-caches that code and
exceeds the limit. Where code nests, the build units must nest too.

A consequence: nested units build serially. Parallelism has to come from units that
genuinely do not share code -- across *components* (atmosphere / ocean / land) and across
*configurations* -- not across phases within one model.

## Type fidelity is what makes the cache usable

Cached specializations are keyed on concrete types. A workload that builds a different type
than the production run yields an image the run cannot use, which looks exactly like
"precompilation didn't help". Two traps found the hard way:

* `FLOAT_TYPE` and `topography` are set in the *coupler* config, not the atmosphere one, so
  a standalone atmosphere workload silently defaults to `Float32` and `Val{:NoWarp}`.
* The gravity-wave model is a type parameter of `AtmosModel`, so a workload that omits
  orographic gravity waves caches specializations nothing uses.

`configs/atmos_warmup.yml` therefore replicates the overrides
`ClimaCoupler.SimCoordinator` applies (`surface_setup`, `albedo_model`, `log_progress`,
`output_dir_style`, `dt`, `toml`), and `CouplerWarmupCommon.check_matches_coupled` fails
loudly on drift.

## The warmup must be the production configuration, briefly

Resolution changes array sizes rather than types, so a cheap low-resolution warmup looks
like it should work. It does not, and three separate floors were found the hard way:

* **Type fidelity** -- model selections must match production exactly (above).
* **Physical validity** -- at `z_elem=10` the coupled state degrades within a few steps
  until `@assert minimum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ 1`
  (`FieldExchanger.jl`) fails at step 3. Confirmed resolution-independent horizontally:
  identical failure at `h_elem` 8 and 16, while production `z_elem=63` ran 5 steps cleanly.
* **Schedule coverage** -- code behind a schedule is only compiled if the schedule fires.
  `orchestrate_diagnostics` skips `compute_field!` entirely when inactive, and the same
  applies to radiation (`dt_rad`), cloud fraction, gravity-wave drag and checkpointing.

So the configuration here is production's, at production resolution, run briefly via `run!`
with every cadence shortened enough to fire. Periods are safe to shorten: `time_to_period`
maps everything except `Nmonths` to `Dates.Millisecond`, so the schedule type is preserved.
ClimaLand accepts only `{1months, 10days, 1days, 1hours, 30mins}`, making `30mins` its floor
and setting `t_end`.

## Building

These packages are developed from this repository rather than registered, so the manifest
is generated locally (and git-ignored, because it records absolute paths):

```bash
julia --project=experiments/warmup -e '
  using Pkg
  Pkg.develop(path=".")
  Pkg.develop([PackageSpec(path=joinpath("experiments","warmup",p)) for p in
    ("CouplerWarmupCommon","AtmosWarmupSetup","AtmosWarmupExplicit","AtmosWarmupStep",
     "OceanWarmup","LandWarmup","SeaIceWarmup","CMIPWarmup")])'

CLIMACOMMS_DEVICE=CUDA JULIA_NUM_PRECOMPILE_TASKS=1 \
  julia --project=experiments/warmup -e 'using Pkg; Pkg.precompile()'

julia --project=experiments/warmup experiments/warmup/check_warmup.jl
```

**Always run `check_warmup.jl` afterwards.** Two failures are invisible otherwise:

* A workload that threw is caught and recorded in `workload_report.txt`, so Pkg sees a
  successful precompile and builds everything above it anyway. Worse, the resulting image
  looks current, so Pkg will not rebuild it later -- the degradation is sticky until the
  image is deleted.
* A link that overflowed the 2 GiB limit happens *after* the workload returns, so no
  workload-level error handling sees it. The symptom is a stale image beside `status = ok`.

A failure low in the DAG makes everything above it re-cache that content, which can push
those packages over the limit too. Fix the lowest failure first, delete the affected images,
and rebuild.

### Worker count is bounded by GPU memory

`JULIA_NUM_PRECOMPILE_TASKS=1` above is deliberate. Each worker builds a full
production-resolution simulation, which at `z_elem=63` holds ~9 GB on the GPU; two
concurrent workers exhaust a 15 GB card and several packages fail with `Out of GPU memory`.
Host memory is not the binding constraint (workers peak near 28 GB against 125 GB). Raise
the count only if the GPU can hold that many simulations at once.

`CLIMACOMMS_DEVICE=CUDA` is required, or the cached specializations are for CPU-backed
fields and a GPU run gains nothing.

## Do not put these in a developer's environment

The workload packages sit *downstream* of ClimaAtmos. Editing ClimaAtmos invalidates them,
and rebuilding is expensive -- one observed cascade ran 3.7 h at 28 GB RSS before being
killed. A developer iterating on, say, microphysics would pay that before every run, making
iteration **worse**, not better.

Precompilation workloads speed up *launches of unedited code*: CI, longruns, production.
They do nothing for someone editing the model, because package images are invalidated at
package granularity. Keep them in this separate environment.

## Cost at run time

Loading the full stack (7 images, ~7.6 GB) takes ~99 s, against 17 s to load the same
packages without it. That is the per-launch cost, set against a setup and first step
measured in tens of minutes. The number of images does not matter -- splitting is free at
run time; total bytes is what costs.

What precompilation cannot reach: **GPU kernel compilation**. CUDA.jl 5.11 cannot persist
device code, so ~2390 kernels are compiled on every launch warm or cold -- 36 % of a cold
first step. Only CUDA.jl 6 changes this, and it will add cubins to images already near the
2 GiB limit, so that migration and this splitting work are coupled.
