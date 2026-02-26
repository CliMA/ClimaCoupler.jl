# Models

ClimaCoupler.jl implements a few simple surface models for use in coupled
simulations. These models are provided in the `Models` module
and provide basic ocean and sea ice representations that can be used when more
sophisticated component models are not needed or available.

## Slab Ocean Model

The slab ocean model (`SlabOceanSimulation`) is a simple energy balance model
for the ocean surface layer. It solves the following energy equation:

```math
(h \cdot \rho \cdot c) \frac{dT}{dt} = -F_{\text{turb\_energy}} + (1 - \alpha) \cdot SW_d + LW_d - LW_u
```

where:
- ``h`` is the depth of the ocean slab [m]
- ``\rho`` is the density of the ocean [kg/m³]
- ``c`` is the specific heat capacity of the ocean [J/(kg·K)]
- ``T`` is the ocean surface temperature [K]
- ``F_{\text{turb\_energy}}`` is the turbulent energy flux (latent + sensible
  heat) [W/m²]
- ``\alpha`` is the ocean albedo
- ``SW_d`` is the downwelling shortwave radiation [W/m²]
- ``LW_d`` is the downwelling longwave radiation [W/m²]
- ``LW_u`` is the upwelling longwave radiation [W/m²] (calculated as
  ``\epsilon \sigma T^4``)

The model assumes a well-mixed surface layer of fixed depth and computes the
temperature evolution based on the net energy flux into the ocean. The ocean
surface temperature is used to compute upwelling longwave radiation and surface
fluxes.

**Usage**: The slab ocean model is used in:
- `SlabplanetMode` simulations
- `SlabplanetAquaMode` simulations

The slab ocean model should not be used with a sea ice model, as it does not
account for ice-ocean interactions.

## Prescribed Ocean Model

The prescribed ocean model (`PrescribedOceanSimulation`) uses observed sea
surface temperature (SST) data that is read from a file at each timestep. The
SST is prescribed and does not evolve based on energy fluxes, making this a
"data-driven" ocean representation.

The model includes:
- Prescribed SST from observational data (e.g., HadISST)
- Ocean roughness parameterization following
  [NOAA-GFDL ice_param](https://github.com/NOAA-GFDL/ice_param/blob/main/ocean_rough.F90#L47)
- Wind-dependent albedo calculation using ClimaAtmos regression functions
- Surface properties (roughness lengths, albedo) that can be updated based on
  atmospheric conditions

Since the SST is prescribed, this model does not solve any prognostic
equations. It serves as a boundary condition for the atmosphere, providing
realistic ocean surface temperatures for atmospheric simulations.

**Usage**: The prescribed ocean model is used in:
- `AMIPMode` simulations
- `SubseasonalMode` simulations

These simulation types are designed for atmospheric model evaluation and
prediction, where realistic SSTs are needed but ocean dynamics are not the
focus.

## Prescribed Sea Ice Model

The prescribed sea ice model (`PrescribedIceSimulation`) uses observed sea ice
concentration (SIC) data while solving an energy balance equation for the ice
surface temperature. The model solves:

```math
(h \cdot \rho \cdot c) \frac{dT_{\text{bulk}}}{dt} = -F_{\text{turb\_energy}} + (1 - \alpha) \cdot SW_d + \epsilon \cdot (LW_d - LW_u) + F_{\text{conductive}}
```

with the conductive flux:

```math
F_{\text{conductive}} = \frac{k_{\text{ice}}}{h} (T_{\text{base}} - T_{\text{sfc}})
```

where:
- ``T_{\text{bulk}}`` is the bulk ice temperature [K] (prognostic variable)
- ``T_{\text{sfc}}`` is the ice surface temperature [K] (diagnostic,
  extrapolated from ``T_{\text{bulk}}`` and ``T_{\text{base}}``)
- ``T_{\text{base}}`` is the prescribed temperature at the ice base (sea water
  temperature) [K]
- ``k_{\text{ice}}`` is the thermal conductivity of sea ice [W/(m·K)]
- ``h`` is the ice thickness [m]
- Other variables are as defined for the slab ocean model

The sea ice concentration is read from observational data (e.g., HadISST) and
updated at each timestep. The ice temperature is constrained to remain at or
below the freezing point of seawater. The model assumes a linear temperature
profile through the ice layer, allowing the surface temperature to be computed
from the bulk temperature and base temperature.

This formulation follows the approach of
[Holloway and Manabe (1971)](https://journals.ametsoc.org/view/journals/mwre/99/5/1520-0493_1971_099_0335_socbag_2_3_co_2.xml),
where sea ice concentration and thickness are prescribed, and the model solves
for the ice temperature.

**Usage**: The prescribed sea ice model is used in:
- `AMIPMode` simulations
- `SubseasonalMode` simulations

These simulation types use prescribed sea ice concentration along with
prescribed ocean SSTs to provide realistic surface boundary conditions for
atmospheric models.


## Models in Extensions

More scientifically-complex component models are provided as package extensions, keeping them out of
ClimaCoupler.jl's base dependencies. The models available in each extension are described below.

### Land models (`ClimaCouplerClimaLandExt`)

This extension is loaded when [ClimaLand.jl](https://github.com/CliMA/ClimaLand.jl) is available
and provides two land model options:

#### `BucketSimulation`

A simple bucket hydrology and energy balance model (`ClimaLand.Bucket.BucketModel`) following the
[Manabe (1969)](https://journals.ametsoc.org/view/journals/mwre/97/11/1520-0493_1969_097_0739_catoc_2_3_co_2.xml)
bucket scheme. This model tracks soil moisture and temperature using a single-layer bucket scheme
and is suitable for lower-cost simulations where a full land model is not needed.

- **Prognostic variables**: temperature, subsurface water content, snow water equivalent
- **Spatial discretization**: finite difference, independent columns on the sphere
- **Time stepping**: explicit additive Runge–Kutta
- **Land surface albedo**: bare ground prescribed from files; snow albedo constant

#### `ClimaLandSimulation`

The full integrated land model (`ClimaLand.LandModel`), which includes soil, canopy,
and snow sub-models. This is the recommended land model for high-fidelity simulations.

- **Prognostic variables**:
  - Soil: internal energy, water content, ice content, carbon content
  - Canopy: temperature, water content
  - Snow: snow water equivalent, snow liquid water, energy
- **Spatial discretization**: finite difference, independent columns on the sphere
- **Time stepping**: mixed implicit/explicit (IMEX) additive Runge–Kutta

For full model equations and further documentation, see the
[ClimaLand.jl documentation](https://clima.github.io/ClimaLand.jl/stable/).

Note that the integrated land model requires spatially-varying parameters that are ill-defined
over ocean/sea ice regions, so it cannot be used in simulations that ignore the land/sea mask
(e.g. terraplanet).

### Ocean and sea ice models (`ClimaCouplerCMIPExt`)

This extension is loaded when [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl),
[ClimaOcean.jl](https://github.com/CliMA/ClimaOcean.jl),
[ClimaSeaIce.jl](https://github.com/CliMA/ClimaSeaIce.jl),
and KernelAbstractions.jl are available. It provides:

- **`OceananigansSimulation`**: A full ocean circulation model implemented in Oceananigans.jl and
  initialized from EN4 or ECCO reanalysis data via ClimaOcean.jl. The ocean evolves
  prognostically in response to atmospheric fluxes. For details on the ocean model physics, see
  the [Oceananigans.jl documentation](https://clima.github.io/OceananigansDocumentation/stable/).

- **`ClimaSeaIceSimulation`**: A thermodynamic sea ice model implemented in ClimaSeaIce.jl. Sea ice
  concentration and thickness evolve prognostically, with ice-ocean heat fluxes computed at the
  ocean-ice interface at each coupling step.

!!! warning "Oceananigans and ClimaSeaIce must be used together"
    `OceananigansSimulation` and `ClimaSeaIceSimulation` are tightly coupled through a shared
    ocean-ice interface with flux exchange and must always be used as a pair. They cannot be combined
    with the simpler ocean or sea ice models provided in `src/`.

### Atmosphere model (`ClimaCouplerClimaAtmosExt`)

This extension is loaded when [ClimaAtmos.jl](https://github.com/CliMA/ClimaAtmos.jl) is
available and provides:

#### `ClimaAtmosSimulation`

The CliMA atmospheric model, which handles radiation, turbulence, convection, microphysics,
and dynamics. This is the only atmosphere model currently supported by ClimaCoupler.jl.

**Dynamical core:**
- Equations: non-hydrostatic and fully compressible
- Prognostic variables: density, velocity components, total energy, total specific humidity
- Spatial discretization: spectral element in the horizontal, finite difference in the vertical, cubed-sphere geometry
- Time stepping: implicit-explicit (IMEX) additive Runge–Kutta

**Parameterizations:**
- Radiation: [RRTMGP](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621)
  (Pincus et al. 2019) — a scheme based on RRTM for GCM applications
- Turbulence and convection: Diagnostic or prognostic Eddy-Diffusivity Mass-Flux (EDMF) scheme with
  prognostic turbulent kinetic energy
  ([Tan et al. 2018](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017MS001162),
  [Cohen et al. 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002162),
  [Lopez-Gomez et al. 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002161))
- Microphysics: 0-moment scheme — cloud condensate removed in-situ with a constant timescale; 1-moment scheme
- Surface fluxes: Monin-Obukhov similarity theory with constant roughness lengths over land
  and ocean

## Model Selection

The choice of component models is determined by the simulation type (`mode_name`).
For a full description of each simulation type and its component model configuration,
see [Available Simulation Types](@ref). Models are selected via configuration files or
command-line arguments; see the [Input](@ref) documentation for details.
