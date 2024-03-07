# ClimaCoupler.jl

*Coupling CliMA Models*

```@meta
CurrentModule = ClimaCoupler
```

ClimaCoupler.jl provides a means to couple [CliMA](https://github.com/CliMA) 
model components. It is designed to provide a flexible way to map boundary fluxes
of quantities, like moisture and heat, that leave one component model
(for example the atmosphere) to boundary fluxes of another component model
(for example the ocean model).
Functionality includes:
- coupled system time stepping control that integrates fluxes in time for sharing
  between components with differing time steps and/or time stepping schemes.
- support for mapping import and export boundary information between components
  so that fluxes of properties transferred between components are conserved.

The ClimaCoupler supports coupling components that are all within the same process
or coupling components (using MPI) that are running on different processes.

| ![Coupler Scheme](images/cplsetup.png) |
|:--:|
| *ClimaCoupler.jl allows for independent development of interchangeable component models.* |


## Current AMIP components

### Atmosphere

Dynamical core:
- Equation: non-hydrostatic and fully compressible
- Prognostic variables: Density, velocity components, total energy, total specific humidity
- Spatial discretization: Spectral element in the horizontal, finite difference in the vertical, cubed sphere
- Time stepping: Implicit-explicit additive Runge–Kutta method

Radiation: A scheme based on RRTM for General circulation model applications—Parallel (RRTMGP)

Convection and turbulence: Diagnostic Eddy-diffusivity Mass-Flux scheme with prognostic turbulent kinetic energy

Microphysics: 0-moment scheme, where cloud condensate is removed with a constant timescale

Surface fluxes: A scheme based on Monin-Obukhov similarity theory, with constant roughness lengths over land and ocean

Orographic gravity wave drag: None

Non-orographic gravity wave drag: None

Aerosols and chemistry: None

### Land

Dynamical core:
- Equation: Following Manabe bucket hydrology scheme
- Prognostic variables: Temperature, water content, snow water content
- Spatial discretization: Finite difference, single column
- Time stepping: Explicit additive Runge–Kutta method

Land surface albedo: 
- Bare ground: Prescribed from files
- Snow: Constant

### Sea ice

Thermodynamics: 0-layer model, with prognostic ice thickness and ice surface temperature

### Coupling

Sequential coupling. Fluxes over a heterogeneous surface are calculated using the averaged surface properties.

```@docs
ClimaCoupler
```
