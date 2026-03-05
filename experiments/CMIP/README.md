# CMIP Experiments

This directory contains the Julia environment and driver scripts for running
CMIP simulations, which require the full dependency stack including Oceananigans,
ClimaOcean, ClimaSeaIce, and KernelAbstractions.

This environment is a superset of `experiments/AMIP/` and can also run AMIP or
slabplanet configurations, though the AMIP environment is preferred for those
since it has fewer dependencies and faster precompilation.

## Running

```bash
julia --project=experiments/CMIP experiments/CMIP/run_simulation.jl --config_file config/ci_configs/cmip_oceananigans_climaseaice.yml
```

## Supported simulation modes

- `cmip` -- fully coupled atmosphere-ocean-land-seaice (Oceananigans + ClimaSeaIce)

For AMIP, slabplanet, and subseasonal modes, use `experiments/AMIP/` instead.
