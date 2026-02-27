# ClimaCoupler.jl

*Coupling CliMA Models*

```@meta
CurrentModule = ClimaCoupler
```

ClimaCoupler.jl provides a flexible infrastructure to connect independent
component models of the Earth system (atmosphere, land, ocean, and sea ice).
This is done by exchanging
information at their boundaries, typically in the form of fluxes.

Key functionality of ClimaCoupler.jl includes:
- **Time stepping**: Coupled system time stepping control allowing for component models to use
  various timesteps individually.
- **Flux calculation**: Flux calculation between the atmosphere and each surface models
  utilizing SurfaceFluxees.jl.
- **Information exchange**: A generic interface to retrieve information from component models that is
  needed to compute surface fluxes, and to update component models with the
  resulting fluxes.
- **Remapping**: Remapping capabilities to support component models running on different
  grid discretizations and with different resolutions.
- **Parallelization**: High-performance computing is enabled via use of GPUs and MPI for distributed
  computing, using ClimaCore.jl to handle under-the-hood numerics and discretizations.

ClimaCoupler.jl is developed as part of the
[Climate Modeling Alliance (CliMA)](https://clima.caltech.edu/) project.

| ![Coupler Scheme](images/cplsetup.png) |
|:--:|
| *ClimaCoupler.jl allows for independent development of interchangeable component models.* |

## Package structure

ClimaCoupler.jl is organized as a set of modules, each with a distinct role in the coupling.

## Coupler interface

### Core modules
- [Interfacer](@ref): Defines the interface for coupling component models, i.e. functions that must
  be extended to use a component model with ClimaCoupler.jl.
- [Input](@ref): Handles parsing command-line inputs and loading configuration files, as well
  as verifying the validity of user-provided options.
- [SimCoordinator](@ref): Coordinates execution of coupled simulations, providing `step!` and
  `run!` to advance the simulation through time.
- [TimeManager](@ref): Handles simulation dates and times, including callback scheduling and
  support for both Float64 and integer time (`ITime`) types.
- [FieldExchanger](@ref): Provides functions required to exchange fields between the atmosphere
  and surface component models at each coupling timestep, building on top of the Interfacer module.
- [FluxCalculator](@ref): Computes turbulent surface fluxes between the atmosphere and each surface
  model using SurfaceFluxes.jl.
- [Checkpointer](@ref): Saves and restores simulation checkpoints, allowing long simulations to be
  restarted as needed.
- [ConservationChecker](@ref): Logs global conservation of energy and water across the coupled
  system, for use with closed-system configurations (i.e. slabplanet simulations).
- [Utilities](@ref): Provides shared helper functions and constants used across modules, including
  device setup, MPI context, and output directory management.

### Provided component models

- [Models](@ref): Implements simple built-in surface models (slab ocean, prescribed ocean,
  prescribed sea ice) for use when more sophisticated component models are not needed.

### Simulation output

- [SimOutput](@ref): Provides utilities for configuring diagnostics output, benchmarking analysis,
  and loading simulation and observational data.
- [Plotting](@ref): Provides visualization tools for simulation output, including diagnostic plots,
  leaderboard comparisons between simulation output and observations, calibration plots, and conservation plots.

### Parameter calibration

- [CalibrationTools](@ref): Provides utilities for setting up model calibration experiments,
  including data loaders for ERA5 observational data.


## Extensions

In addition to the modules located in `src/`, ClimaCoupler.jl contains package extensions
located in `ext/`. These primarily extend the coupling interface for particular component models,
and also provide plotting utilities. Placing these pieces in extensions allows the top-level
package environment of ClimaCoupler.jl to remain streamlined, while also supporting more
complex setups.

The available extensions are:

- `ClimaCouplerClimaAtmosExt`: Extends the coupler interface for
  [ClimaAtmos.jl](https://github.com/CliMA/ClimaAtmos.jl), implementing field exchanges, flux
  computation, and checkpointing for the CliMA atmosphere model.
- `ClimaCouplerClimaLandExt`: Extends the coupler interface for
  [ClimaLand.jl](https://github.com/CliMA/ClimaLand.jl), supporting both the bucket land model
  and the integrated land model.
- `ClimaCouplerCMIPExt`: Extends the coupler interface for CMIP-class coupled simulations,
  integrating [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) as the ocean model
  and [ClimaSeaIce.jl](https://github.com/CliMA/ClimaSeaIce.jl) as the sea ice model.
- `ClimaCouplerMakieExt`: Implements all [Plotting](@ref) functions using Makie.jl,
  CairoMakie.jl, ClimaCoreMakie.jl, and GeoMakie.jl. Loaded automatically when these packages
  are available.
- `ClimaCouplerCMIPMakieExt`: Extends the Makie plotting support to handle
  Oceananigans.jl fields when Oceananigans is used as the ocean component model.

Details about the component model extensions can be found in the section [Models in Extensions](@ref),
and details about the plotting extensions can be found in [Plotting](@ref).

## Documentation structure

- **[Running a simulation](@ref)**: How to launch and configure a simulation from the
  command line, REPL, or as an installed package.
- **[Available component models](@ref)**: Describes the built-in surface models (slab
  ocean, prescribed ocean, prescribed sea ice) and the more complex models provided by
  package extensions.
- **[Available simulation types](@ref)**: Explains each supported `mode_name` option
  (AMIP, CMIP, slabplanet variants, subseasonal) and the component models each uses.
- **Coupler interface**: This section includes one page per module in `src/`, each documenting
  that module's public API and functions.
- **Simulation output**: Documentation for saving diagnostics, visualization, leaderboard
  comparisons against observations, and performance benchmarking.
- **Examples**: Walkthroughs of running specific experiments.
- **[CalibrationTools](@ref)**: Utilities for setting up parameter calibration experiments
  using ERA5 observational data.
- **Developer docs**: Guidelines for contributing to ClimaCoupler.jl, and practical tips
  for debugging package incompatibilities, numerical instabilities, and software errors.
