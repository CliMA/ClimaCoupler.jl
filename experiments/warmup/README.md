# Precompilation workloads

Packages whose `@compile_workload` runs real model code, so that the inference and codegen
for setup and the first step land in a package image instead of being redone on every launch.

**These are deliberately *not* dependencies of `experiments/CMIP` or `experiments/AMIP`.**
See "Do not put these in a developer's environment" below.

## Why

Nothing in the CliMA stack ships `@compile_workload`s, so each launch re-infers and re-emits
the whole setup and step path on one core. Measured on a coupled CMIP configuration
(prognostic EDMF, Float64, CUDA):

| phase | wall | of which compilation |
|---|---:|---:|
| package loading | 17 s | — |
| `CoupledSimulation(...)` | 61 min | 84 % |
| first `step!` | 19–23 min | 91 % |
| second `step!` | 0.2 s | 0 % |

Startup is ~99 % single-core: Julia serializes codegen process-wide, so the only parallelism
available is across *processes*, which is what `Pkg.precompile()` provides.

These figures, and the image sizes below, were measured on Julia 1.11 / CUDA 5.11 /
ClimaCore 0.16. They have not been re-measured since the move to Julia 1.13 and CUDA 6.

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

Image sizes at production resolution: 1260 / 1205 / 925 / 1240 / 1332 / 108 / 1554 MiB.
`SeaIceWarmup` is small because `OceanWarmup` is beneath it and it caches only the
incremental sea-ice code.

The config loading, `WARMUP_H_ELEM` override and `workload_report.txt` writer every package
uses live in `warmup_utils.jl`, which each package `include`s. It is a file, not a package,
so sharing it adds no edge to the graph above.

## The dependency edges are required, not an optimisation

A package image cannot exceed ~2 GiB. Above that, code in the image references the
serialized blob with 32-bit ELF relocations (`R_X86_64_PC32`, `GOTPC32`) that cannot encode
the displacement. Over the limit `lld` errors out; just under it, the link succeeds silently
and the image segfaults on load inside `jl_restore_package_image_from_stream`.

The edges must mirror the **call graph**. `step!` calls `remaining_tendency!`, so
`AtmosWarmupStep` depends on `AtmosWarmupExplicit`; without that edge it re-caches the
tendency code and exceeds the limit. The same applies to `SeaIceWarmup` and `OceanWarmup`.
Where code nests, the build units must nest too.

A single unsplit workload produces one 3.05 GiB object that `lld` refuses to link.

A consequence: nested units build serially. Parallelism has to come from units that
genuinely do not share code -- across *components* (atmosphere / ocean / land) and across
*configurations* -- not across phases within one model.

## Type fidelity is what makes the cache usable

Cached specializations are keyed on concrete types. A workload that builds a different type
than the production run yields an image the run cannot use, which looks exactly like
"precompilation didn't help". Two settings are easy to get wrong:

* `FLOAT_TYPE` and `topography` are set in the *coupler* config, not the atmosphere one, so
  a standalone atmosphere workload defaults to `Float32` and `Val{:NoWarp}`.
* The gravity-wave model is a type parameter of `AtmosModel`, so a workload that omits
  orographic gravity waves caches specializations nothing uses.

`configs/atmos_warmup.yml` therefore replicates the overrides
`ClimaCoupler.SimCoordinator` applies (`surface_setup`, `albedo_model`, `log_progress`,
`output_dir_style`, `dt`, `toml`), and `CouplerWarmupCommon.check_matches_coupled` fails
loudly on drift.

## The warmup must be the production configuration, briefly

A cheap low-resolution warmup looks like it should work. Three constraints rule it out:

* **Vertical resolution is part of the type.** A field's data layout is
  `ClimaCore.DataLayouts.VIJFH{FT, Nv, ...}`, where `Nv` is the number of levels, so an
  image built at one `z_elem` is not usable by a run at another. Measured against an image
  built at `z_elem = 79`, `ClimaAtmos.get_simulation` compiled 0.2 s at 79, 250.6 s at 78
  and 243.8 s at 63 -- one level of difference costs essentially the whole cache.
  Horizontal resolution is *not* affected: `h_elem` changes the number of elements, which
  is a runtime dimension.
* **Type fidelity** -- model selections must match production exactly (above).
* **Physical validity** -- the warmup has to be a configuration that actually runs. Well
  below production vertical resolution the coupled state degrades within a few steps until
  `@assert minimum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ 1`
  (`FieldExchanger.jl`) fails.
* **Schedule coverage** -- code behind a schedule is only compiled if the schedule fires.
  `orchestrate_diagnostics` skips `compute_field!` entirely when inactive, and the same
  applies to radiation (`dt_rad`), cloud fraction, gravity-wave drag and checkpointing.

### Building on a small GPU

The coupled workloads hold the atmosphere, ocean and sea ice on the device at once, which
at production resolution exceeds a 15 GB card. `WARMUP_H_ELEM` lowers the *horizontal*
resolution, which cuts memory without costing cache validity, since `h_elem` is not part of
the type. Images built at `h_elem = 8` serve a production run at 16.

```bash
CLIMACOMMS_DEVICE=CUDA JULIA_NUM_PRECOMPILE_TASKS=1 WARMUP_H_ELEM=8 \
  julia --project=experiments/warmup -e 'using Pkg; Pkg.precompile()'
```

There is deliberately no equivalent for `z_elem`: it is the `Nv` type parameter, so
lowering it produces images no production run can use. The value is recorded in
`workload_report.txt`.

So the configuration here is production's, at production resolution, run briefly via `run!`
with every cadence shortened enough to fire. Periods are safe to shorten: `time_to_period`
maps everything except `Nmonths` to `Dates.Millisecond`, so the schedule type is preserved.
ClimaLand accepts only `{1months, 10days, 1days, 1hours, 30mins}`, making `30mins` its floor
and setting `t_end`.

## Building

The packages are developed from this repository rather than registered. Their paths live in
`Project.toml`'s `[sources]` block and are relative, so the manifest is committed and the
environment instantiates without a separate `Pkg.develop` step:

```bash
CLIMACOMMS_DEVICE=CUDA JULIA_NUM_PRECOMPILE_TASKS=1 \
  julia --project=experiments/warmup -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

```bash
julia --project=experiments/warmup experiments/warmup/check_warmup.jl
```

### Versions must match `experiments/CMIP`

A cached specialization is keyed on the package versions it was built from, so the warmup
environment has to resolve to exactly what the run environment resolves to. The manifest
here is therefore seeded from `experiments/CMIP/Manifest-v1.13.toml` and only the
warmup-only extras are resolved on top of it. Left to float, it picks up newer releases
(observed: ClimaCore 1.0.8 against CMIP's 1.0.5) and the cache silently stops applying.

Re-seed after any change to the CMIP manifest:

```bash
cp experiments/CMIP/Manifest-v1.13.toml experiments/warmup/Manifest-v1.13.toml
julia --project=experiments/warmup -e 'using Pkg; Pkg.resolve()'
```

**Always run `check_warmup.jl` afterwards.** Two failures are invisible otherwise:

* A workload that threw is caught and recorded in `workload_report.txt`, so Pkg sees a
  successful precompile and builds everything above it anyway. Worse, the resulting image
  looks current, so Pkg will not rebuild it later -- the degradation is sticky until the
  image is deleted.
* A link that overflowed the 2 GiB limit happens *after* the workload returns, so no
  workload-level error handling sees it. The symptom is a stale image beside `status = ok`.

A failure low in the DAG makes everything above it re-cache that content, which can push
those packages over the limit too. Fix the lowest failure first, delete the affected images
(a caught failure leaves a current-looking image that Pkg will not rebuild), and rebuild.

### Worker count is bounded by GPU memory

`JULIA_NUM_PRECOMPILE_TASKS=1` above is deliberate. Each worker builds a full
production-resolution simulation holding ~9 GB on the GPU, so two concurrent workers exhaust
a 15 GB card. Host memory is not the binding constraint (workers peak near 28 GB). Raise the
count only if the GPU can hold that many simulations at once.

`CLIMACOMMS_DEVICE=CUDA` is required, or the cached specializations are for CPU-backed
fields and a GPU run gains nothing.

## Do not put these in a developer's environment

The workload packages sit *downstream* of ClimaAtmos and ClimaCore, so editing either
invalidates them. A full rebuild is ~95 min, against ~2 min to rebuild the 22 packages a
ClimaCore change invalidates without them -- so a batch `Pkg.precompile` and run workflow on
model source pays roughly +95 min per edit. Breakeven is about two coupled runs per edit.

Precompilation workloads speed up *launches of unedited code*: CI, longruns, production,
and calibration (where each persistent worker pays startup once). Iterating on model source
in a REPL with Revise should be roughly neutral, since Revise patches methods in the live
session rather than rebuilding images, but that has not been measured.

Keeping them in this separate environment means nobody pays for them who is not using them.

## Cost at run time

Loading the full stack (7 images, ~7.6 GB) takes ~99 s, against 17 s to load the same
packages without it. That is the per-launch cost, set against a setup and first step
measured in tens of minutes. The number of images does not matter -- splitting is free at
run time; total bytes is what costs.

**GPU kernel compilation** used to be out of reach: CUDA.jl 5.11 could not persist device
code, so ~2390 kernels were compiled on every launch warm or cold -- 36 % of a cold first
step. On the current stack (Julia 1.13, CUDA 6.4.1, GPUCompiler 2.10) it can be:
`GPUCompiler.can_persist_results` is true for a CUDA job, because `supports_relocatable_ir()`
holds from Julia 1.13 and `CUDACore` requests `relocation_lowering = :patch`.

The consequence for this directory is that persisted device code *adds* cubins to images
that already peaked at 1554 MiB against the ~2 GiB ceiling below. Image sizes need
re-measuring, and the split may need to go further.
