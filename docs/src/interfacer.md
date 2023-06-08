# Interfacer

This module contains functions for defining the interface for coupling component models, as well as stub objects that contain prescribed fields.

## Coupled Simulation
- `CoupledSimulation` (`cs`, currently in Utilities - TODO) stores info for ESM run. We require that each `cs` contains four (`atmos_sim`, `land_sim`, `ocean_sim` and `ice_sim`) components. While this requirement will not be eventually needed, for the time being, if a simulation surface type is not needed for a given run, it should be initialized with `SurfaceStub` with a zero `area_fracion`. The `atmos_sim` should always be specified.

## Component model simulations
- all Simulations that are not the `CoupledSimulation` fall under `ComponentModelSimulation`
- the current version requires that there is:
    - one `AtmosModelSimulation`
    - one or more `SurfaceModelSimulation`s, which require the following adapter functions:
        ```
        get_field(sim::SurfaceModelSimulation, ::Val{:area_fraction}) = ...
        get_field(sim::SurfaceModelSimulation, ::Val{:surface_temperature}) = ...
        get_field(sim::SurfaceModelSimulation, ::Val{:albedo}) = ...
        get_field(sim::SurfaceModelSimulation, ::Val{:roughness_momentum}) = ...
        get_field(sim::SurfaceModelSimulation, ::Val{:roughness_buoyancy}) = ...
        get_field(sim::SurfaceModelSimulation, ::Val{:beta}) = ...
        update_field!(sim::SurfaceModelSimulation, ::Val{:area_fraction}, field::Fields.Field) = ...
        update_field!(sim::SurfaceModelSimulation, ::Val{:surface_temperature}, field::Fields.Field) = ...
        ```
        - these adapter functions, to be defined in the component models' init files (preferably in their own repositories), allow the coupler to operate without having to assume particular data structures of the underlying component models. This allows easy swapping of model components, as well as a stable source code with coupler-specific unit tests.

## Prescribed conditions
- `SurfaceStub` is a `SurfaceModelSimulation`, but it only contains required data in `<surface_stub>.cache`, e.g., for the calculation of surface fluxes through a prescribed surface state.  The above adapter functions are already predefined for `SurfaceStub`, with the cache variables specified as:
```
get_field(sim::SurfaceStub, ::Val{:area_fraction}) = sim.cache.area_fraction
get_field(sim::SurfaceStub, ::Val{:surface_temperature}) = sim.cache.T_sfc
get_field(sim::SurfaceStub, ::Val{:albedo}) = sim.cache.α
get_field(sim::SurfaceStub, ::Val{:roughness_momentum}) = sim.cache.z0m
get_field(sim::SurfaceStub, ::Val{:roughness_buoyancy}) = sim.cache.z0b
get_field(sim::SurfaceStub, ::Val{:beta}) = sim.cache.beta
```
with the corresponding `update_field!` functions
```
function update_field!(sim::SurfaceStub, ::Val{:area_fraction}, field::Fields.Field)
    sim.cache.area_fraction .= field
end
function update_field!(sim::SurfaceStub, ::Val{:surface_temperature}, field::Fields.Field)
    sim.cache.T_sfc .= field
end
```

## Interfacer API

```@docs
    ClimaCoupler.Interfacer.ComponentModelSimulation
    ClimaCoupler.Interfacer.AtmosModelSimulation
    ClimaCoupler.Interfacer.SurfaceModelSimulation
    ClimaCoupler.Interfacer.SurfaceStub
    ClimaCoupler.Interfacer.name
    ClimaCoupler.Interfacer.get_field
    ClimaCoupler.Interfacer.update_field!
```
