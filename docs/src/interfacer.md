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
- constructor: construct and return an instance of the `ComponentModelSimulation`,
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

### ComponentModelSimulation - optional functions
- `get_model_prog_state(::ComponentModelSimulation)`: A function that
returns the state vector of the simulation at its current state. This
is used for checkpointing the simulation.

- `get_field(::ComponentModelSimulation, ::Val{property})`:
Default `get_field` functions are provided for `energy` and `water` fields,
described in the table below.
These quantities are used to track conservation, and the defaults
return `nothing`. To check conservation throughout a simulation, these
functions must be extended for all models being run.

| Coupler name | Description                                   | Units  | Default value |
|--------------+-----------------------------------------------+--------+---------------|
| `energy`     | vertically integrated energy per surface area | J m⁻²  | `nothing`     |
| `water`      | vertically integrated water per surface area  | kg m⁻² | `nothing`     |

- `add_coupler_fields!(coupler_field_names, ::ComponentModelSimulation)`:
A set of default coupler exchange fields is initialized for each coupled simulation,
but depending on the component models being run, additional coupler fields may
be required. For example, the integrated land model requires the concentration
of atmospheric CO2 for photosynthesis calculations, but the slab ocean does not.
`add_coupler_fields!` is extended for any component model simulation that requires
coupler fields in addition to the defaults, allowing us to allocate space for and
exchange the extra fields only when necessary.
  - Any additional fields specified here will likely also require an `update_field!`
method for this component model to pass the value from the coupler to the component.

The default coupler exchange fields are the following, defined in
`default_coupler_fields()` in the Interfacer module:

| Coupler name      | Description                                         | Units      |
|-------------------+-----------------------------------------------------+------------|
| `z0m_sfc`         | momentum roughness length                           | m          |
| `z0b_sfc`         | buoyancy roughness length                           | m          |
| `beta`            | factor to scale evaporation from the surface        | -          |
| `emissivity`      | surface emissivity                                  | -          |
| `F_turb_energy`   | turbulent energy flux                               | W m⁻²      |
| `F_turb_moisture` | turbulent moisture flux                             | kg m⁻² s⁻¹ |
| `F_turb_ρτxz`     | turbulent momentum flux in the zonal direction      | kg m⁻¹ s⁻² |
| `F_turb_ρτyz`     | turbulent momentum flux in the meridional direction | kg m⁻¹ s⁻² |
| `P_liq`           | liquid precipitation                                | kg m⁻² s⁻¹ |
| `P_snow`          | snow precipitation                                  | kg m⁻² s⁻¹ |
| `temp1`           | a surface field used for intermediate calculations  | -          |
| `temp2`           | a surface field used for intermediate calculations  | -          |

- `update_sim!(::ComponentModelSimulation, csf)`: A
function to update each of the fields of the component model simulation
that are updated by the coupler. ClimaCoupler.jl provides defaults of
this function for both `AtmosModelSimulation` and
`SurfaceModelSimulation` that update each of the fields expected by
the coupler. This function will need to be extended for any model
that requires additional fields (specified via `add_coupler_fields!`).

### AtmosModelSimulation - required functions
In addition to the functions required for a general
`ComponentModelSimulation`, an `AtmosModelSimulation` requires the
following functions to retrieve and update its fields.
- `get_field(::AtmosModelSimulation. ::Val{property})`: This getter
function returns the value of the field property for the simulation
in its current state. For an `AtmosModelSimulation`, it must be extended
for the following properties:

| Coupler name                | Description                                                               | Units      |
|-----------------------------+---------------------------------------------------------------------------+------------|
| `air_density`               | air density of the atmosphere                                             | kg m⁻³     |
| `height_int`                | height at the first internal model level                                  | m          |
| `height_sfc`                | height at the surface                                                     | m          |
| `liquid_precipitation`      | liquid precipitation at the surface                                       | kg m⁻² s⁻¹ |
| `radiative_energy_flux_sfc` | net radiative flux at the surface                                         | W m⁻²      |
| `radiative_energy_flux_toa` | net radiative flux at the top of the atmosphere                           | W m⁻²      |
| `snow_precipitation`        | snow precipitation at the surface                                         | kg m⁻² s⁻¹ |
| `turbulent_energy_flux`     | aerodynamic turbulent surface fluxes of energy (sensible and latent heat) | W m⁻²      |
| `turbulent_moisture_flux`   | aerodynamic turbulent surface fluxes of energy (evaporation)              | kg m⁻² s⁻¹ |
| `thermo_state_int`          | thermodynamic state at the first internal model level                     |            |
| `uv_int`                    | horizontal wind velocity vector at the first internal model level         | m s⁻¹      |

- `update_field!(::AtmosModelSimulation. ::Val{property}, field)`:
A function to update the value of property in the component model
simulation, using the values in `field`. This update should
be done in place. If this function isn't extended for a property,
that property will remain constant throughout the simulation
and a warning will be raised.
This function is expected to be extended for the
following properties, and may also be extended for any additional
properties needed by a component model.

| Coupler name             | Description                                              | Units |
|--------------------------+----------------------------------------------------------+-------|
| `emissivity`             | surface emissivity                                       |       |
| `surface_direct_albedo`  | bulk direct surface albedo over the whole surface space  |       |
| `surface_diffuse_albedo` | bulk diffuse surface albedo over the whole surface space |       |
| `surface_temperature`    | temperature over the combined surface space              | K     |
| `turbulent_fluxes`       | turbulent fluxes                                         | W m⁻² |

- `calculate_surface_air_density(atmos_sim::Interfacer.AtmosModelSimulation, T_sfc::ClimaCore.Fields.Field)`:
A function to return the air density of the atmosphere simulation
extrapolated to the surface, with units of [kg m⁻³].

ClimaAtmos should also add the following coupler fields for Monin-Obukhov similarity theory:
| Coupler name    | Description       | Units  |
|-----------------+-------------------+--------|
| `ustar`         | friction velocity | m s⁻¹  |
| `L_MO`          | Obukhov length    | m      |
| `buoyancy_flux` | flux of buoyancy  | m⁻²s⁻³ |


### AtmosModelSimulation - required functions to run with the ClimaLandSimulation

Coupling with the integrated `ClimaLandSimulation` requires the following functions, in addition
to the functions required for coupling with a general `SurfaceModelSimulation`.

- `get_field(::AtmosModelSimulation. ::Val{property})`:
This getter function must be extended
for the following properties:

| Coupler name        | Description                                                    | Units   |
|---------------------+----------------------------------------------------------------+---------|
| `air_pressure`      | air pressure at the bottom cell centers of the atmosphere      | Pa      |
| `air_temperature`   | air temperature at the bottom cell centers of the atmosphere   | K       |
| `cos_zenith`        | cosine of the zenith angle                                     |         |
| `co2`               | global mean co2                                                | ppm     |
| `diffuse_fraction`  | fraction of downwards shortwave flux that is direct            |         |
| `specific_humidity` | specific humidity at the bottom cell centers of the atmosphere | kg kg⁻¹ |
| `LW_d`              | downwards longwave flux                                        | W m⁻²   |
| `SW_d`              | downwards shortwave flux                                       | W m⁻²   |

Note that `air_temperature`, `air_pressure`, `cos_zenith`, `co2`, `diffuse_fraction`, `LW_d` and
`SW_d` will not be present in a `ClimaAtmosSimulation` if the model is setup with no radiation.
Because of this, a `ClimaAtmosSimulation` must have radiation if running with the full `ClimaLand` model.

### SurfaceModelSimulation - required functions
Analogously to the `AtmosModelSimulation`, a `SurfaceModelSimulation`
requires additional functions to those required for a general `ComponentModelSimulation`.
- `get_field(::SurfaceModelSimulation, ::Val{property})`: This getter
function returns the value of the field property for the simulation at
the current time. For a `SurfaceModelSimulation`, it must be extended
for the following properties:

| Coupler name             | Description                                                    | Units   |
|--------------------------+----------------------------------------------------------------+---------|
| `area_fraction`          | fraction of the simulation grid surface area this model covers |         |
| `roughness_buoyancy`     | aerodynamic roughness length for buoyancy                      | m       |
| `roughness_momentum`     | aerodynamic roughness length for momentum                      | m       |
| `surface_direct albedo`  | bulk direct surface albedo                                     |         |
| `surface_diffuse albedo` | bulk diffuse surface albedo                                    |         |
| `surface_humidity`       | surface humidity                                               | kg kg⁻¹ |
| `surface_temperature`    | surface temperature                                            | K       |


- `update_field!(::SurfaceModelSimulation, ::Val{property}, field)`:
A function to update the value of property in the component model
simulation, using the values in `field` passed from the coupler
This update should be done in place. If this function
isn't extended for a property,
that property will remain constant throughout the simulation
and a warning will be raised.
This function is expected to be extended for the
following properties, and may also be extended for any additional
properties needed by a component model.

| Coupler name                                  | Description                                                                  | Units      |
|-----------------------------------------------+------------------------------------------------------------------------------+------------|
| `air_density`                                 | surface air density                                                          | kg m⁻³     |
| `area_fraction`                               | fraction of the simulation grid surface area this model covers               |            |
| `liquid_precipitation`                        | liquid precipitation at the surface                                          | kg m⁻² s⁻¹ |
| `radiative_energy_flux_sfc` OR `LW_d`, `SW_d` | net radiative flux at the surface OR downward longwave, shortwave radiation  | W m⁻²      |
| `snow_precipitation`                          | snow precipitation at the surface                                            | kg m⁻² s⁻¹ |
| `turbulent_energy_flux`                       | aerodynamic turbulent surface fluxes of energy (sensible and latent heat)    | W m⁻²      |
| `turbulent_moisture_flux`                     | aerodynamic turbulent surface fluxes of energy (evaporation)                 | kg m⁻² s⁻¹ |


### SurfaceModelSimulation - optional functions
- `get_field(::SurfaceModelSimulation, ::Val{property})`:
For some quantities, default `get_field` functions are provided, which may be
overwritten or used as-is. These currently include the following:

| Coupler name  | Description                                                               | Units | Default value |
|---------------+---------------------------------------------------------------------------+-------+---------------|
| `beta`        | factor that scales evaporation based on its estimated level of saturation |       |             1 |
| `emissivity`  | measure of how much energy a surface radiates                             |       |             1 |
| `height_disp` | displacement height relative to the surface                               | m     |             0 |


- `update_turbulent_fluxes!(::ComponentModelSimulation, fields::NamedTuple)`:
This function updates the turbulent fluxes of the component model simulation
at this point in horizontal space. The values are updated using the energy
and moisture turbulent fluxes stored in fields which are calculated by the
coupler.

### Prescribed surface conditions - SurfaceStub
- `SurfaceStub` is a `SurfaceModelSimulation`, but it only contains
required data in `<surface_stub>.cache`, e.g., for the calculation
of surface fluxes through a prescribed surface state. The above
adapter functions are already predefined for `AbstractSurfaceStub`,
which is extended by `SurfaceStub`
in the `surface_stub.jl` file, with
the cache variables specified as:
```
get_field(sim::AbstractSurfaceStub, ::Val{:area_fraction}) = sim.cache.area_fraction
get_field(sim::AbstractSurfaceStub, ::Val{:beta}) = sim.cache.beta
get_field(sim::AbstractSurfaceStub, ::Val{:roughness_buoyancy}) = sim.cache.z0b
get_field(sim::AbstractSurfaceStub, ::Val{:roughness_momentum}) = sim.cache.z0m
get_field(sim::AbstractSurfaceStub, ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}}) = sim.cache.α
get_field(sim::AbstractSurfaceStub, ::Val{:surface_humidity}) = TD.q_vap_saturation_generic.(sim.cache.thermo_params, sim.cache.T_sfc, sim.cache.ρ_sfc, sim.cache.phase)
get_field(sim::AbstractSurfaceStub, ::Val{:surface_temperature}) = sim.cache.T_sfc
```
and with the corresponding `update_field!` functions
```
function update_field!(sim::AbstractSurfaceStub, ::Val{:air_density}, field)
    sim.cache.ρ_sfc .= field
end
function update_field!(sim::AbstractSurfaceStub, ::Val{:area_fraction}, field::ClimaCore.Fields.Field)
    sim.cache.area_fraction .= field
end
function update_field!(sim::AbstractSurfaceStub, ::Val{:surface_temperature}, field::ClimaCore.Fields.Field)
    sim.cache.T_sfc .= field
end
```

## Interfacer API
```@docs
    ClimaCoupler.Interfacer.CoupledSimulation
    ClimaCoupler.Interfacer.AtmosModelSimulation
    ClimaCoupler.Interfacer.SurfaceModelSimulation
    ClimaCoupler.Interfacer.ComponentModelSimulation
    ClimaCoupler.Interfacer.AbstractSurfaceStub
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
    ClimaCoupler.Interfacer.AbstractSimulation
    ClimaCoupler.Interfacer.AbstractSimulationMode
```
