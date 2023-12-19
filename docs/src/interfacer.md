# Interfacer

This module contains functions for defining the interface for coupling component models, as well as stub objects that contain prescribed fields.

## Coupled Simulation
- `CoupledSimulation` (`cs`) stores info for ESM run. We require that each `cs` contains four (`atmos_sim`, `land_sim`, `ocean_sim` and `ice_sim`) components. While this requirement will not be eventually needed, for the time being, if a simulation surface type is not needed for a given run, it should be initialized with `SurfaceStub` with a zero `area_fracion`. The `atmos_sim` should always be specified.

## Component model simulations
- all Simulations that are not the `CoupledSimulation` fall under `ComponentModelSimulation`
- the current version requires that there is:
    - one `AtmosModelSimulation`
    - one or more `SurfaceModelSimulation`s, which require the following adapter functions:


## `AtmosModelSimulation` requirements

## `SurfaceModelSimulation` requirements
### "Getter" functions (`get_field`)
- Each `SurfaceModelSimulation` is required to implement the following "getter"
functions, which retrieve a value from the model and provide it to the coupler.
This requires accessing the model's internal structure, so they should be
implemented in an `init` function in the model's own repository to abstract
model specifics away from the coupler.
    ```
    get_field(sim::SurfaceModelSimulation, ::Val{:surface_temperature}) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:surface_humidity}) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:roughness_momentum}) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:roughness_buoyancy}) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:beta}) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:albedo}) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:area_fraction}, field::Fields.Field) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:air_density}, field::Fields.Field) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:energy}, field::Fields.Field) = ...
    get_field(sim::SurfaceModelSimulation, ::Val{:water}, field::Fields.Field) = ...
    ```
    - these adapter functions, to be defined in the component models' init files (preferably in their own repositories), allow the coupler to operate without having to assume particular data structures of the underlying component models. This allows easy swapping of model components, as well as a stable source code with coupler-specific unit tests.
    - note that these functions are currently defined in component models' init
    files in `experiments/AMIP`, but these will be moved to the components'
    own repositories eventually.

### Update functions (`update_field!`)
- Each `SurfaceModelSimulation` is also required to implement the following
update methods. These should overwrite the model's value for a particular
quantity with the value stored in `field`, again accessing the internals of
the particular model:
    ```
    update_field!(sim::SurfaceModelSimulation, ::Val{:turbulent_energy_flux}, field) = ...
    update_field!(sim::SurfaceModelSimulation, ::Val{:turbulent_moisture_flux}, field) = ...
    update_field!(sim::SurfaceModelSimulation, ::Val{:radiative_energy_flux}, field) = ...
    update_field!(sim::SurfaceModelSimulation, ::Val{:liquid_precipitation}, field) = ...
    update_field!(sim::SurfaceModelSimulation, ::Val{:snow_precipitation}, field)
    update_field!(sim::SurfaceModelSimulation, ::Val{:air_density}, field)
    ```

### Timestepping functions
- Each component model is required to provide methods for the following
functions so that the coupler is able to timestep each component model:
    ```
    step!(sim::SurfaceModelSimulation, t) = ...
    reinit!(sim::SurfaceModelSimulation) = ...
    ```

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
