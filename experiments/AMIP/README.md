# AMIP Experiments

This directory contains the Julia environment and driver scripts for running
AMIP, slabplanet, and subseasonal simulations. It does **not** include
dependencies needed for CMIP simulations (Oceananigans, ClimaSeaIce,
KernelAbstractions, ConservativeRegridding). For CMIP, see `experiments/CMIP/`.

## Running

Run these commands from the repository root with Julia 1.10 or 1.11, the versions
this environment's manifests are pinned to. Install the environment once:

```bash
julia --project=experiments/AMIP -e 'using Pkg; Pkg.instantiate()'
```

Then run a simulation:

```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl
```

To run with a specific configuration file:

```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl --config_file config/ci_configs/slabplanet_default.yml
```

## Supported simulation modes

- `amip` -- prescribed SST/SIC atmosphere benchmark
- `slabplanet` / `slabplanet_aqua` / `slabplanet_terra` -- idealized configurations
- `subseasonal` -- short-range forecasts from ERA5 initial conditions

See the top-level [README.md](https://github.com/CliMA/ClimaCoupler.jl/blob/main/README.md) and
`experiments/CMIP/README.md` for the CMIP configuration.
