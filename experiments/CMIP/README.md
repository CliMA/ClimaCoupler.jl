# CMIP Experiments

This directory contains the Julia environment and driver scripts for running
CMIP simulations, which require the full dependency stack including Oceananigans,
ClimaSeaIce, KernelAbstractions, ConservativeRegridding, and Adapt.

This environment is a superset of `experiments/AMIP/` and can also run AMIP or
slabplanet configurations, though the AMIP environment is preferred for those
since it has fewer dependencies and faster precompilation.

## Running

Run these commands from the repository root with Julia 1.10 or 1.11, the versions
this environment's manifests are pinned to. Install the environment once:

```bash
julia --project=experiments/CMIP -e 'using Pkg; Pkg.instantiate()'
```

Then run a simulation:

```bash
julia --project=experiments/CMIP experiments/CMIP/run_simulation.jl --config_file config/ci_configs/cmip_oceananigans_climaseaice.yml
```

## Supported simulation modes

- `cmip` -- fully coupled atmosphere-ocean-land-seaice (Oceananigans + ClimaSeaIce)

For AMIP, slabplanet, and subseasonal modes, use `experiments/AMIP/` instead.
