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

## Model Selection

The choice of which models to use is typically determined by the simulation
mode:

- **AMIP/Subseasonal modes**: Use `prescribed` ocean and `prescribed` sea ice
- **CMIP mode**: Uses more sophisticated models (`oceananigans` ocean,
  `clima_seaice` sea ice)
- **Slabplanet modes**: Use `slab` ocean and no sea ice (or `nothing` for both
  ocean and ice in terra mode)

The models are selected via configuration files or command-line arguments when
setting up a `CoupledSimulation`. For more information about how to select these
component models, please see the `Input` module documentation.
