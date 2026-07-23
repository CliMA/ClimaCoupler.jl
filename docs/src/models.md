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
- Ocean COARE3 roughness parameterization, which accounts for the
  dynamic response of ocean surface roughness to wind speed and wave state
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

This model requires the `ClimaCouplerClimaAtmosExt` extension to be loaded in order to use
the wind-dependent albedo calculation, which relies on the `surface_albedo_direct` and `surface_albedo_diffuse`
functions from the `ClimaAtmos` module.
If the extension is not loaded, users must define their own `set_albedos!` function for the `PrescribedOceanSimulation`.

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

## Eisenman-Zhang Sea Ice Model

The Eisenman-Zhang sea ice model (`EisenmanIceSimulation`) is a thermodynamic
0-layer sea ice model over a slab ocean mixed layer, based on the
[Semtner (1976)](https://doi.org/10.1175/1520-0485(1976)006<0379:AMFTTG>2.0.CO;2)
model and later refined by
[Eisenman & Wettlaufer (2009)](https://doi.org/10.1073/pnas.0806887106) and
[Zhang et al. (2021)](https://doi.org/10.1029/2021MS002671)
(whose prior implementation can be found on
[GitHub](https://github.com/sally-xiyue/fms-idealized/blob/sea_ice_v1.0/exp/sea_ice/srcmods/mixed_layer.f90)).
Unlike the prescribed sea ice model, the ice thickness evolves prognostically,
with ice growing or melting in response to the surface energy imbalance, and
the surface can transition between ice-covered and ice-free states. The model
assumes no snow coverage.

The prognostic variables are:
- ice thickness ``h_i``
- ocean mixed layer temperature ``T_{ml}``
- surface temperature ``T_s``
(plus an accumulated basal energy used for energy bookkeeping).
Ice cover is a binary mask: the surface is ice-covered wherever ``h_i > 0``.

### Ice covered

In ice-covered conditions the **ice thickness** ``h_i`` evolves as

```math
L_i \frac{dh_i}{dt} = F_{atm} - F_{base} - Q
```

where:
- ``L_i`` is the volumetric latent heat of fusion of ice (default: ``3.0e8`` J m^-3)
- ``t`` is model time in seconds
- ``F_{atm}`` is the total energy flux from the surface to the atmosphere (positive upwards)
- ``F_{base}`` is the basal heat flux from the ocean mixed layer into the ice
- ``Q`` is the prescribed ocean q-flux (positive = oceanic heat source; see
  [Ocean q-flux](@ref) below)

The energy flux into the atmosphere can be expanded as:
```math
F_{atm} = F_{rad} + F_{turb\_energy}
```

and further where:

```math
F_{rad} = \epsilon (\sigma T_s^4 - LW_d) - (1 - \alpha) SW_d
```

whereas the basal heat flux is taken as:
```math
F_{base} = C_0(T_{ml} - T_{base})
```

where:
- ``C_0`` is the linear basal heat transfer coefficient (default: ``120`` W m^-2 K^-1)
- ``T_{base}`` is the ice base temperature (default: ``273.16`` K, equal to the
  freezing temperature ``T_{freeze}``)

The **mixed layer temperature** under ice responds to the basal flux and the
ocean q-flux (the mixed layer only exchanges heat with the atmosphere when
ice-free):

```math
\rho_w c_w h_{ml} \frac{dT_{ml}}{dt} = -(F_{base} - Q)
```

The **ice surface temperature** ``T_s`` is obtained by balancing the total surface energy flux ``F_{atm}(T_s)`` against the conductive heat flux through the ice slab,
``F_{ice} = k_i (T_{base} - T_s) / h_i`` with ``k_i = 2`` W m^-1 K^-1 (default):

```math
F_{atm} = F_{ice} = k_i \frac{T_{base} - T_s}{h_i}
```

Solving using one Newton iteration (sufficient at the current spatial and temporal
resolution — see Semtner, 1976):

```math
T_s^{t+1} = T_s^{t} + \frac{- F_{atm}^t + k_i (T_{base} - T_s^t) / h_i^{t+1}}{k_i/h_i^{t+1} + \partial F_{atm}^t / \partial T_s^t}
```

where the conductive flux is evaluated with the updated ice thickness
``h_i^{t+1}`` and the current surface temperature ``T_s^t``. The updated
``T_s`` is capped at the freezing point (the ice surface stores no energy). The
derivative ``\partial F_{atm} / \partial T_s = 4 \epsilon \sigma T_s^3``
contains only the radiative contribution: the turbulent flux derivative
``\partial F_{\text{turb}} / \partial T_s`` is no longer provided by the
coupler (its finite-difference machinery was removed in
[#1284](https://github.com/CliMA/ClimaCoupler.jl/pull/1284)), so the
turbulent flux is treated explicitly in the Newton update.

### Ice free

In ice-free conditions, the mixed layer assumes the standard slab representation

```math
\rho_w c_w h_{ml} \frac{dT_{ml}}{dt} = -(F_{atm} - Q)
```

and ``T_s = T_{ml}``. The default mixed layer parameters are depth
``h_{ml} = 1`` m, density ``\rho_w = 1020`` kg m^-3, and heat capacity
``c_w = 4000`` J kg^-1 K^-1.

**Frazil ice formation**: the mixed layer is not allowed to cool below
``T_{freeze}``; the energy deficit is instead used to grow ice.

**Transition to ice-free conditions**: if the updated ``h_i`` would be
negative, the ice thickness is set to zero and the surplus energy warms the
mixed layer.

### Ocean q-flux

The forcing term ``Q`` (`ocean_qflux` in the model cache) represents a
prescribed oceanic heat convergence into the column, with the convention that
positive ``Q`` is a heat source. In ice-free conditions it enters the
mixed-layer budget alongside ``F_{atm}`` (warming the mixed layer); in
ice-covered conditions it both warms the mixed layer (offsetting ``F_{base}``)
and thins the ice, as in the equations above. The accumulated q-flux
``Q \, t`` is also included in the column energy bookkeeping returned by
`Interfacer.get_field(sim, Val(:energy))`.

Note that ``Q`` is currently always zero in coupled runs: it is initialized to
zero in `eisenman_state_init` and no `Interfacer.update_field!` method or
coupler exchange pathway sets it, so a nonzero q-flux can only be imposed by
writing to `sim.integrator.p.ocean_qflux` directly (as the unit tests do).

**Usage**: select with `ice_model: "eisenman"`. Since the model carries its
own mixed layer, it should not be paired with a separate ocean model
overlapping the same surface area.

### Potential extensions

- add an `ice_area_fraction` adjustment (e.g., assuming a minimal thickness of
  sea ice, below which the grid area becomes part ice and part ocean)


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

## Julia Environments

ClimaCoupler provides two Julia environments under `experiments/`, each targeting
a different set of simulation modes:

### `experiments/AMIP`

Use this environment for `amip`, `slabplanet`, `slabplanet_aqua`, `slabplanet_terra`,
and `subseasonal` modes. It omits the Oceananigans-related dependencies (`Oceananigans`,
`ClimaOcean`, `ClimaSeaIce`, `KernelAbstractions`), resulting in a lighter environment
with faster precompilation.

```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl
```

### `experiments/CMIP`

Use this environment for `cmip` mode, which requires a prognostic Oceananigans ocean
model. It is a superset of the AMIP environment with additional ocean and sea-ice
dependencies.

```bash
julia --project=experiments/CMIP experiments/CMIP/run_simulation.jl --config_file config/ci_configs/cmip_oceananigans_climaseaice.yml
```

See the [Input](@ref) documentation for the full list of `mode_name` options.
