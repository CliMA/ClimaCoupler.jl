# SimOutput

The `SimOutput` module provides utilities for output operations including
diagnostics and performance analysis.

## Diagnostics Setup

The diagnostics setup function configures default diagnostics for AMIP simulations,
which uses ClimaDiagnostics.jl to save variables throughout the course of a simulation.

For more information about diagnostics in ClimaCoupler, including how to customize which
variables to save, how often, and with which reductions, see the [Diagnostics](@ref) documentation.

### Functions

```@docs
SimOutput.diagnostics_setup
```

## Benchmarking Analysis

The benchmark analysis functions help compare performance metrics (e.g. SYPD, or simulated years per day)
between different simulation runs. This information is formatted into a table format using
PrettyTables.jl, and is sent to Slack automatically each time the "benchmarks" buildkite pipeline is run.

### Functions

```@docs
SimOutput.get_benchmark_args
SimOutput.get_run_info
SimOutput.append_table_data
```

## Simulation and Observation Data

These functions provide dictionaries mapping diagnostic variable names to preprocessed simulation and observational
data for use in leaderboard comparisons. They are used by [`compute_leaderboard`](@ref) to load
and prepare variables for comparison against observations.

For more information about adding variables to the leaderboard, see the [Leaderboard](@ref) documentation.

### Functions

```@docs
SimOutput.get_sim_var_dict
SimOutput.get_obs_var_dict
```
