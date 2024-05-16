# AMIP Log

## Current AMIP components

### Atmosphere

Dynamical core:
- Equation: non-hydrostatic and fully compressible
- Prognostic variables: Density, velocity components, total energy, total specific humidity
- Spatial discretization: Spectral element in the horizontal, finite difference in the vertical, cubed sphere
- Time stepping: Implicit-explicit additive Runge–Kutta method

Parameterizations:

- Radiation: A scheme based on RRTM for General circulation model applications—Parallel (RRTMGP) (Pincus et al. 2019)

- Convection and turbulence: Diagnostic Eddy-diffusivity Mass-Flux (EDMF) scheme with prognostic turbulent kinetic energy.
EDMF is a unified parameterization for turbulence and convection (Tan et al. 2018, Cohen et al. 2020, Lopez-Gomez et al. 2020). The grid is decomposed into convective updrafts and the turbulent environment.
Updraft properties are calculated from mass, momentum and energy conservation of the updraft.
Mass and momentum exchange between the updrafts and the environment, and
turbulent mixing in the environment, are represented with physical closures. Currently, only one updraft is used.

- Microphysics: 0-moment scheme, where cloud condensate is removed in-situ with a constant timescale

- Surface fluxes: A scheme based on Monin-Obukhov similarity theory, with constant roughness lengths over land and ocean

- Orographic gravity wave drag: None

- Non-orographic gravity wave drag: None

- Aerosols and chemistry: None

### Land

Dynamical core:
- Equation: Following Manabe bucket hydrology scheme (Manabe 1969)
- Prognostic variables: Temperature, water content, snow water content
- Spatial discretization: Finite difference, multiple independent columns on the sphere
- Time stepping: Explicit additive Runge–Kutta method

Land surface albedo:
- Bare ground: Prescribed from files
- Snow: Constant

### Sea ice

Thermodynamics: 0-layer model, with prognostic ice surface temperature and fixed ice thickness

### Coupling

Sequential coupling. Fluxes over a heterogeneous surface are calculated using the averaged surface properties.

## References

[Manabe 1969](https://journals.ametsoc.org/view/journals/mwre/97/11/1520-0493_1969_097_0739_catoc_2_3_co_2.xml): Climate and ocean circulation. I. The atmospheric circulation and the hydrology of the earth's surface.

[Pincus et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621): Balancing Accuracy, Efficiency, and Flexibility in Radiation Calculations for Dynamical Models

[Tan et al. 2018](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017MS001162): An Extended Eddy-Diffusivity Mass-Flux Scheme for Unified Representation of Subgrid-Scale Turbulence and Convection

[Cohen et al. 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002162): Unified Entrainment and Detrainment Closures for Extended Eddy-Diffusivity Mass-Flux Schemes

[Lopez-Gomez et al. 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002161): A Generalized Mixing Length Closure for Eddy-Diffusivity Mass-Flux Schemes of Turbulence and Convection