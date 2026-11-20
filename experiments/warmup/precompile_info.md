# Precompilation: measurements and constraints

Hard numbers behind the design of `experiments/warmup` and `experiments/CMIP/startup`. The
code and its docstrings state the rules; this file records the evidence for them. Each
section names the stack it was measured on, because most of these numbers move with Julia,
CUDA and package versions.

## Where startup time goes

Coupled CMIP configuration (prognostic EDMF, Float64, CUDA), no warmup images.
Julia 1.11 / CUDA 5.11 / ClimaCore 0.16.

| phase | wall | of which compilation |
|---|---:|---:|
| package loading | 17 s | — |
| `CoupledSimulation(...)` | 61 min | 84 % |
| first `step!` | 19–23 min | 91 % |
| second `step!` | 0.2 s | 0 % |

* Startup is ~99 % single-core: Julia serializes codegen process-wide, so the only
  parallelism available is across processes, which is what `Pkg.precompile()` provides.
* Within the 61-minute setup: inference 18.8 min, LLVM codegen 34.5 min. SnoopCompile sees
  inference only, so it misses the larger half; `phase_probe.jl` samples C frames instead.
* On a run this size `Profile.print` adds ~1 h of post-processing, and sampling correlated
  with hangs when two heavy jobs ran concurrently. Hence `cmip_startup_nosample.jl`.

## Package image size limit (~2 GiB)

Code in a package image references the serialized blob through 32-bit ELF relocations
(`R_X86_64_PC32`, `GOTPC32`). Above ~2 GiB the displacement cannot be encoded: `lld` refuses
to link, and just under the limit the link can succeed and the image segfault on load in
`jl_restore_package_image_from_stream`.

* A single unsplit workload produced one 3.05 GiB object that `lld` refused to link.
* Hence the split into packages, with dependency edges that follow the call graph so each
  package reuses, rather than re-caches, the code beneath it.

## Image sizes

| package | Julia 1.11 / CUDA 5.11 | `check_warmup.jl` floor (Julia 1.13 / CUDA 6.4) |
|---|---:|---:|
| `AtmosWarmupSetup` | 1260 MiB | 300 MiB |
| `AtmosWarmupExplicit` | 1205 MiB | 150 MiB |
| `AtmosWarmupStep` | 925 MiB | 150 MiB |
| `OceanWarmup` | 1240 MiB | 200 MiB |
| `LandWarmup` | 1332 MiB | 200 MiB |
| `SeaIceWarmup` | 108 MiB | 25 MiB |
| `CMIPWarmup` | 1554 MiB | 75 MiB |

* The whole stack shrank about 2.7× between Julia 1.11 / CUDA 5.11 and Julia 1.13 / CUDA 6.4,
  even though the latter also persists GPU kernels into the image (~89 KiB each).
* On Julia 1.13: `SeaIceWarmup` ~80 MiB, `CMIPWarmup` ~219 MiB. Both are small because the
  packages beneath them carry the components. A run against that stack reached its first
  step with 5.7 s of setup compilation.
* The floors are set to catch an image that cached nothing, not to track ordinary drift.

## PrecompileTools ≥ 1.3 on Julia ≥ 1.12

PrecompileTools 1.2 caches no externally-owned code on Julia ≥ 1.12, while the workload
still runs and reports `status = ok`.

* Same workload: 0.02 MiB cached with 1.2, 143.55 MiB with 1.3.
* Observed in a full build: `AtmosWarmupExplicit` at 2 MiB where ~1200 MiB was expected.

## Load order: what invalidates the cache

A package loaded after a workload's code is compiled invalidates that code, and GPU kernels
attached to the invalidated `CodeInstance`s are discarded with it. The warmup packages
therefore import everything `experiments/CMIP/code_loading.jl` loads, even where nothing
calls it.

* Makie stack absent from the workload: `ClimaAtmos.get_simulation` compiled 0.2 s in a
  session without Makie and 227.1 s in one with it; coupled setup went 4.4 s → 1153.1 s.
* Coupled packages absent from the atmosphere workload: step 1 linked its 260 stepping
  kernels in 0.9 s with only the atmosphere chain loaded, and recompiled all 260 at a cost
  of 128.9 s once the coupled packages were loaded.

## Resolution and types

A field's data layout is `ClimaCore.DataLayouts.VIJFH{S, Nv, Nij, Nij, ...}`: it carries
the vertical level count `Nv` and quadrature size, not the number of horizontal elements.

* `h_elem` is type-neutral. Against an image built at `h_elem = 16`, `get_simulation` at
  `h_elem = 8` compiled 0.2 s, the same as at 16.
* `z_elem` is not. Against an image built at `z_elem = 79`, `get_simulation` compiled 0.2 s
  at 79, 250.6 s at 78 and 243.8 s at 63.

## Configuration and version drift

Several model selections are type parameters, so a warmup config that drifts from production
yields an image the run cannot use, and the run silently recompiles.

* `FLOAT_TYPE` and `topography` are set in the coupler config, not the atmosphere one, so a
  standalone atmosphere workload defaulted to `Float32` and `Val{:NoWarp}`. Caught by
  `CouplerWarmupCommon.check_matches_coupled`.
* Production switched `orographic_gravity_wave` from `raw_topo` to `raw_topo_online` on
  2026-10-02; the warmup config kept the old value for a month. Now caught by
  `check_config_drift.jl`.
* A fresh resolve under Julia 1.13 floated to ClimaCore 1.0.8, ClimaAtmos 0.42.18 and
  ClimaLand 1.15.0 against CMIP's 1.0.5 / 0.42.17 / 1.14.1. Hence the manifest is seeded from
  CMIP's, and `check_version_drift.jl` checks it stays so.

## Build and run cost

* Full rebuild of the warmup stack: ~95 min. Without the warmup packages, a ClimaCore edit
  rebuilds 22 packages in ~2 min. Breakeven is about two coupled runs per edit.
* Each build worker holds a production-resolution simulation (~9 GB) on the GPU, so two
  concurrent workers exceed a 15 GB card. Host memory is not the constraint (workers peak
  near 28 GB).
* Loading the full stack (7 images, ~7.6 GB, Julia 1.11) took ~99 s, against 17 s without
  it. Splitting into more images is free at run time; total bytes is what costs.

## GPU kernel persistence

* CUDA.jl 5.11 could not persist device code: ~2390 kernels were compiled on every launch,
  36 % of a cold first step.
* On Julia 1.13 / CUDA 6.4.1 / GPUCompiler 2.10, `GPUCompiler.can_persist_results` is true
  for a CUDA job (`supports_relocatable_ir()` holds from Julia 1.13 and `CUDACore` requests
  `relocation_lowering = :patch`), so kernels are cached in the image.
