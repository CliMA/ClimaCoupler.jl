# Postprocessor

The `Postprocessor` module provides utilities for postprocessing simulation results,
including plotting and performance analysis.

## Diagnostics Setup

The diagnostics setup function configures default diagnostics for AMIP simulations.
For more information about diagnostics in ClimaCoupler, see the [Diagnostics](@ref) documentation.

### Functions

```@docs
Postprocessor.diagnostics_setup
```

## Benchmarking Analysis

The benchmark analysis functions help compare performance metrics (like SYPD - Simulated Years Per Day) between different simulation runs.

### Functions

```@docs
Postprocessor.get_benchmark_args
Postprocessor.get_run_info
Postprocessor.append_table_data
```
