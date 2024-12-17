# Interfacer
This module contains functions that define the interface for coupling
component models, as well as stub objects that contain prescribed fields.
Here we explain each type of component model, and the functions that must
be implemented to use a component model with ClimaCoupler.jl

## Coupled simulation
A `CoupledSimulation` stores info for ESM run and contains each of the
component model simulations. We currently require that each `CoupledSimulation`
contains four components: `atmos_sim`, `land_sim`, `ocean_sim` and `ice_sim`.
If a simulation surface type is not needed for a given run, it should be
initialized with `SurfaceStub` with a zero `area_fracion`.
The `atmos_sim` should always be specified.


## Component simulations
Individual component model simulations fall under `ComponentModelSimulation`,
which together combine to make the `CoupledSimulation`.
We have two types of `ComponentModelSimulations`: `AtmosModelSimulation` and
`SurfaceModelSimulation`. The two have different requirements,
which are detailed below. `SurfaceModelSimulation` is further divided into
`SeaIceModelSimulation`, `LandModelSimulation`, and `OceanModelSimulation`,
representing the 3 currently-supported options for surface models.

### ComponentModelSimulation - required functions
A component model simulation should be implemented as a struct that is a concrete subtype
of a `ComponentModelSimulation`. This struct should contain all of the
information needed to run that simulation.

Each `ComponentModelSimulation` must extend the following functions to be able
to use our coupler. For some existing models, these are defined within
ClimaCoupler.jl in that model’s file in `experiments/ClimaEarth/components/`, but it is preferable
for these to be defined in a model’s own repository. Note that the dispatch
`::ComponentModelSimulation` in the function definitions given below should
be replaced with the particular component model extending these functions.
- `init`: construct and return an instance of the `ComponentModelSimulation`,
and perform all initialization. This function should return a simulation that
is ready to be stepped in the coupled simulation. The interface for this
function varies across component models.
- `name(::ComponentModelSimulation)`: return a string containing the name of
this `ComponentModelSimulation`, which is used for printing information about
component models and writing to checkpoint files.
- `step!(::ComponentModelSimulation, t)`: A function to update the
simulation in-place with values calculate for time `t`. For the
models we currently have implemented, this is a simple wrapper around
the `step!` function implemented in SciMLBase.jl.
- `reinit!(::ComponentModelSimulation)`: A function to restart a simulation
after solving of the simulation has been paused or interrupted. Like
`step!`, this is currently a simple wrapper around the `reinit!` function
of SciMLBase.jl.
- `get_model_prog_state(::ComponentModelSimulation)`: A function that
returns the state vector of the simulation at its current state. This
is used for checkpointing the simulation.

### ComponentModelSimulation - optional functions
- `update_sim!(::ComponentModelSimulation, csf, turbulent_fluxes)`: A
function to update each of the fields of the component model simulation
that are updated by the coupler. ClimaCoupler.jl provides defaults of
this function for both `AtmosModelSimulation` and
`SurfaceModelSimulation` that update each of the fields expected by
the coupler. This function can optionally be extended to include
additional field updates as desired.

### AtmosModelSimulation - required functions
In addition to the functions required for a general
`ComponentModelSimulation`, an `AtmosModelSimulation` requires the
following functions to retrieve and update its fields.
- `get_field(::AtmosModelSimulation. ::Val{property})`: This getter
function returns the value of the field property for the simulation
in its current state. For an `AtmosModelSimulation`, it must be extended
for the following properties:

| Coupler name      | Description | Units |
|-------------------|-------------|-------|
| `air_density`       | air density of the atmosphere | kg m^-3 |
| `air_temperature`   | air temperature at the surface of the atmosphere | K |
| `energy`            | globally integrated energy | J |
| `height_int`        | height at the first internal model level | m |
| `height_sfc`        | height at the surface (only required when using `PartitionedStateFluxes`) | m |
| `liquid_precipitation` | liquid precipitation at the surface | kg m^-2 s^-1 |
| `radiative_energy_flux_sfc` | net radiative flux at the surface | W m^-2 |
| `radiative_energy_flux_toa` | net radiative flux at the top of the atmosphere | W m^-2 |
| `snow_precipitation` | snow precipitation at the surface | kg m^-2 s^-1 |
| `turbulent_energy_flux` | aerodynamic turbulent surface fluxes of energy (sensible and latent heat) | W m^-2 |
| `turbulent_moisture_flux` | aerodynamic turbulent surface fluxes of energy (evaporation) | kg m^-2 s^-1 |
| `thermo_state_int`  | thermodynamic state at the first internal model level | |
| `uv_int`            | horizontal wind velocity vector at the first internal model level | m s^-1 |
| `water`             | globally integrated water | kg |



- `update_field!(::AtmosModelSimulation. ::Val{property}, field)`:
A function to update the value of property in the component model
simulation, using the values in `field`. This update should
be done in place. If this function isn't extended for a property,
that property will remain constant throughout the simulation
and a warning will be raised.
This function is expected to be extended for the
following properties:

| Coupler name      | Description | Units |
|-------------------|-------------|-------|
| `co2`              | global mean co2 | ppm |
| `surface_direct_albedo`   | bulk direct surface albedo over the whole surface space | |
| `surface_diffuse_albedo`   | bulk diffuse surface albedo over the whole surface space | |
| `surface_temperature` | temperature over the combined surface space | K |
| `turbulent_fluxes` | turbulent fluxes (note: only required when using `PartitionedStateFluxes` option - see our `FluxCalculator` module docs for more information) | W m^-2 |

- `calculate_surface_air_density(atmos_sim::Interfacer.AtmosModelSimulation, T_S::ClimaCore.Fields.Field)`:
A function to return the air density of the atmosphere simulation
extrapolated to the surface, with units of [kg m^-3].

### SurfaceModelSimulation - required functions
Analogously to the `AtmosModelSimulation`, a `SurfaceModelSimulation`
requires additional functions to those required for a general `ComponentModelSimulation`.
- `get_field(::SurfaceModelSimulation. ::Val{property})`: This getter
function returns the value of the field property for the simulation at
the current time. For a `SurfaceModelSimulation`, it must be extended
for the following properties:

| Coupler name      | Description | Units |
|-------------------|-------------|-------|
| `air_density`       | surface air density | kg m^-3 |
| `area_fraction`     | fraction of the simulation grid surface area this model covers | |
| `beta`              | factor that scales evaporation based on its estimated level of saturation | |
| `roughness_buoyancy` | aerodynamic roughness length for buoyancy | m |
| `roughness_momentum` | aerodynamic roughness length for momentum | m |
| `surface_direct albedo`    | bulk direct surface albedo | |
| `surface_diffuse albedo`    | bulk diffuse surface albedo | |
| `surface_humidity`  | surface humidity | kg kg^-1 |
| `surface_temperature` | surface temperature | K |

- `update_field!(::SurfaceModelSimulation. ::Val{property}, field)`:
A function to update the value of property in the component model
simulation, using the values in `field` passed from the coupler
This update should be done in place. If this function
isn't extended for a property,
that property will remain constant throughout the simulation
and a warning will be raised.
This function is expected to be extended for the
following properties:

| Coupler name      | Description | Units |
|-------------------|-------------|-------|
| `air_density`       | surface air density | kg m^-3 |
| `area_fraction`     | fraction of the simulation grid surface area this model covers | |
| `liquid_precipitation` | liquid precipitation at the surface | kg m^-2 s^-1 |
| `radiative_energy_flux_sfc` | net radiative flux at the surface | W m^-2 |
| `snow_precipitation` | snow precipitation at the surface | kg m^-2 s^-1 |
| `turbulent_energy_flux` | aerodynamic turbulent surface fluxes of energy (sensible and latent heat) | W m^-2 |
| `turbulent_moisture_flux` | aerodynamic turbulent surface fluxes of energy (evaporation) | kg m^-2 s^-1 |
| `surface_direct_albedo`    | bulk direct surface albedo; needed if calculated externally of the surface model (e.g. ocean albedo from the atmospheric state) | |
| `surface_diffuse_albedo`    | bulk diffuse surface albedo; needed if calculated externally of the surface model (e.g. ocean albedo from the atmospheric state) | |

### SurfaceModelSimulation - optional functions
- `update_turbulent_fluxes!(::ComponentModelSimulation, fields::NamedTuple)`:
This function updates the turbulent fluxes of the component model simulation
at this point in horizontal space. The values are updated using the energy
and moisture turbulent fluxes stored in fields which are calculated by the
coupler. Note that this function is only required when using the
`PartitionedStateFluxes` option of ClimaCoupler.jl. See our `FluxCalculator`
module docs for more information.

### Prescribed surface conditions - SurfaceStub
- `SurfaceStub` is a `SurfaceModelSimulation`, but it only contains
required data in `<surface_stub>.cache`, e.g., for the calculation
of surface fluxes through a prescribed surface state. The above
adapter functions are already predefined for `SurfaceStub`
in the `surface_stub.jl` file, with
the cache variables specified as:
```
get_field(sim::SurfaceStub, ::Val{:air_density}) = sim.cache.ρ_sfc
get_field(sim::SurfaceStub, ::Val{:area_fraction}) = sim.cache.area_fraction
get_field(sim::SurfaceStub, ::Val{:beta}) = sim.cache.beta
get_field(sim::SurfaceStub, ::Val{:energy}) = nothing
get_field(sim::SurfaceStub, ::Val{:roughness_buoyancy}) = sim.cache.z0b
get_field(sim::SurfaceStub, ::Val{:roughness_momentum}) = sim.cache.z0m
get_field(sim::SurfaceStub, ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}}) = sim.cache.α
get_field(sim::SurfaceStub, ::Val{:surface_humidity}) = TD.q_vap_saturation_generic.(sim.cache.thermo_params, sim.cache.T_sfc, sim.cache.ρ_sfc, sim.cache.phase)
get_field(sim::SurfaceStub, ::Val{:surface_temperature}) = sim.cache.T_sfc
get_field(sim::SurfaceStub, ::Val{:water}) = nothing
```
and with the corresponding `update_field!` functions
```
function update_field!(sim::SurfaceStub, ::Val{:air_density}, field)
    sim.cache.ρ_sfc .= field
end
function update_field!(sim::SurfaceStub, ::Val{:area_fraction}, field::ClimaCore.Fields.Field)
    sim.cache.area_fraction .= field
end
function update_field!(sim::SurfaceStub, ::Val{:surface_temperature}, field::ClimaCore.Fields.Field)
    sim.cache.T_sfc .= field
end
```

## Interfacer API
```@docs
    ClimaCoupler.Interfacer.CoupledSimulation
    ClimaCoupler.Interfacer.AtmosModelSimulation
    ClimaCoupler.Interfacer.SurfaceModelSimulation
    ClimaCoupler.Interfacer.ComponentModelSimulation
    ClimaCoupler.Interfacer.SurfaceStub
    ClimaCoupler.Interfacer.float_type
    ClimaCoupler.Interfacer.name
    ClimaCoupler.Interfacer.get_field
    ClimaCoupler.Interfacer.update_field!
    ClimaCoupler.Interfacer.AbstractSlabplanetSimulationMode
    ClimaCoupler.Interfacer.AMIPMode
    ClimaCoupler.Interfacer.SlabplanetMode
    ClimaCoupler.Interfacer.SlabplanetAquaMode
    ClimaCoupler.Interfacer.SlabplanetTerraMode
    ClimaCoupler.Interfacer.SlabplanetEisenmanMode
```

## Interfacer Internal Functions and Types

```@docs
    ClimaCoupler.Interfacer.stub_init
    ClimaCoupler.Interfacer.AbstractSimulation
    ClimaCoupler.Interfacer.AbstractSlabplanetSimulationMode
```
