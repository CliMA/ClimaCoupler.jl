# Performance Tips

## CPU vs GPU performance

CliMA models are primarily tuned for best performance on GPUs, but all simulation types
can run on CPU. The right choice depends on the simulation type:

- **Slabplanet and AMIP** are feasible on CPU and are commonly used for development and
  testing on a laptop or single workstation node.
- **CMIP** is prohibitively slow on CPU due to the combined cost of all the component
  models. While it is technically possible to run on CPU, we strongly recommend using a GPU.

The compute device is controlled by the `device` configuration option. See [Input](@ref)
for details.

### Reference performance numbers

The numbers below were last updated **6 March 2026** and they come from the
buildkite longrun and regular CI pipelines. The turbulence convection (TC)
scheme used in the atmosphere is noted as it has a strong impact on performance.
All runs use 0-moment microphysics, which also greatly affects performance.

| Configuration | Hardware                       | Resolution                              | SYPD  | TC scheme       | Source                      |
|---------------|--------------------------------|-----------------------------------------|-------|-----------------|-----------------------------|
| Slabplanet    | 1 CPU (1 core, Intel Ice Lake) | ~2.8°                                   | 0.034 | ED only         | longrun                     |
| AMIP          | 1 CPU (1 core, Intel Ice Lake) | ~3.8°                                   | 0.276 | vert. diff only | Float64 + hourly checkpoint |
| AMIP          | 1 GPU (A100)                   | ~1.4°                                   | 1.996 | ED only         | longrun                     |
| CMIP          | 1 GPU (A100)                   | ~1.4° (atmos/land), ~1° (ocean/sea ice) | 1.017 | ED only         | longrun                     |

!!! note
    SYPD (simulated years per day) is the standard metric for climate model throughput.
    Higher is better. Resolution is given as approximate grid spacing at the equator.

## Compiler optimization flags

Because `CoupledSimulation` is a complex data structure that contains many component
models and their associated caches, constructing and initializing it can take significant
time. Much of this is spent on Julia's just-ahead-of-time compilation — for background
on how Julia precompilation works, see the
[Julia precompilation tutorial](https://julialang.org/blog/2021/01/precompile_tutorial/).

Julia's `-O` flag controls the compiler optimization level (default: `2`). It accepts
values `0`, `1`, `2`, and `3`, where higher values allow more aggressive specialization
during compilation.

**For very short runs (≤15 coupling steps)**, passing `-O0` is recommended. At this
scale the overall runtime is dominated by compilation time, and reducing the optimization
level significantly shortens the startup cost. This is useful when iterating quickly
on initialization logic or debugging:

```bash
julia -O0 --project=experiments/ClimaEarth experiments/ClimaEarth/run_amip.jl
```

**For longer runs**, the default optimization level (`-O2`) is recommended, as the
compilation overhead becomes negligible relative to the simulation runtime and the
generated code runs faster.
