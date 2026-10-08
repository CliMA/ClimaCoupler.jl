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
| `AtmosWarmupSetup` | `ClimaAtmos.get_simulation` | — |
| `AtmosWarmupExplicit` | `remaining_tendency!` | Setup |
| `AtmosWarmupStep` | implicit solve, Jacobian, timestepper | Setup, Explicit |
| `CouplerWarmupCommon` | shared prologue for component workloads | ClimaCoupler |

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
loudly on drift. Resolution (`h_elem`, `z_elem`) is deliberately small: it changes array
sizes, not types.

## Building

```bash
CLIMACOMMS_DEVICE=CUDA julia --project=experiments/warmup -e 'using Pkg; Pkg.precompile()'
```

`CLIMACOMMS_DEVICE=CUDA` is required, or the cached specializations are for CPU-backed
fields and a GPU run gains nothing. Cap `JULIA_NUM_PRECOMPILE_TASKS` by available memory:
workers peaked near 28 GB, and two concurrent ones nearly exhausted a 125 GB machine.

Each package writes `workload_report.txt` recording whether its workload succeeded, so a
caught failure is visible rather than silently producing an image with nothing useful in it.

## Do not put these in a developer's environment

The workload packages sit *downstream* of ClimaAtmos. Editing ClimaAtmos invalidates them,
and rebuilding is expensive -- one observed cascade ran 3.7 h at 28 GB RSS before being
killed. A developer iterating on, say, microphysics would pay that before every run, making
iteration **worse**, not better.

Precompilation workloads speed up *launches of unedited code*: CI, longruns, production.
They do nothing for someone editing the model, because package images are invalidated at
package granularity. Keep them in this separate environment.

## Cost at run time

Loading the images is cheap: 2466 MiB of package images cost **+7.3 s** over a bare
`import ClimaAtmos` (22.0 s -> 29.3 s), against the ~29 minutes of setup and first step they
replace. The number of images does not matter; splitting is effectively free at run time.
