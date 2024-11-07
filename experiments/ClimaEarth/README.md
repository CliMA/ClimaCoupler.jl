# ClimaEarth Experiments

The main script in the ClimaEarth directory is [run_amip.jl](run_amip.jl). This script can be used to run the AMIP configuration, or a simplified slabplanet configuration, both of which use sequential coupling. Each configuration is explained here.

## AMIP components
AMIP is currently the most complex configuration of the ClimaEarth model.
It runs a ClimaAtmos.jl atmosphere model, ClimaLand.jl bucket land model,
a prescribed ocean model, and a simple thermal sea ice model
Please find more details about each model below.

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

Sea ice concentration is prescribed by the `SIC` input file and changes with time.

### Ocean
The ocean model has no prognostic variables. It stores sea surface temperature, which is
read from the `SST` input file and changes over time.

### References: Atmosphere formulation

[Manabe 1969](https://journals.ametsoc.org/view/journals/mwre/97/11/1520-0493_1969_097_0739_catoc_2_3_co_2.xml): Climate and ocean circulation. I. The atmospheric circulation and the hydrology of the earth's surface.

[Pincus et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621): Balancing Accuracy, Efficiency, and Flexibility in Radiation Calculations for Dynamical Models

[Tan et al. 2018](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017MS001162): An Extended Eddy-Diffusivity Mass-Flux Scheme for Unified Representation of Subgrid-Scale Turbulence and Convection

[Cohen et al. 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002162): Unified Entrainment and Detrainment Closures for Extended Eddy-Diffusivity Mass-Flux Schemes

[Lopez-Gomez et al. 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002161): A Generalized Mixing Length Closure for Eddy-Diffusivity Mass-Flux Schemes of Turbulence and Convection

## Slabplanet components
The Slabplanet configuration is more idealized than the AMIP configuration, but provides valuable insight about conservation and individual model behavior.

### Atmosphere
The atmosphere model used in the slabplanet configuration is the same as the one used for AMIP;
please see the previous section for details.

### Land
The land model used in the slabplanet configuration is the same as the one used for AMIP;
please see the previous section for details.

### Ocean
Thermodynamics: thermal slab model; prognostic sea surface temperature without depth

### Sea ice
The slabplanet configuration does not include sea ice. Instead, ocean is evaluated in areas
that would be covered by sea ice.

### Slabplanet aqua
This configuration is similar to the general "Slabplanet" configuration, except that the
only surface model is the ocean, which is evaluated over the entire surface. There are no land or sea ice models.

### Slabplanet terra
This configuration is similar to the general "Slabplanet" configuration, except that the
only surface model is the land, which is evaluated over the entire surface. There are no ocean or sea ice models.

### Slabplanet Eisenman
This configuration is similar to the general "Slabplanet" configuration, except that the ocean model
is included in the Eisenman sea ice model.

#### Eisenman Sea Ice and Ocean
Thermodynamics: 0-layer model, based on the Semtner 1976 model and later refined by
Eisenman & Wettlaufer (2009) and Zhang et al. (2021).

Prognostic variables: ice height (`h_i`), ocean mixed layer depth (`T_ml`) and surface air temperature (`T_s`).

Note that Eisenman sea ice assumes gray radiation, no snow coverage, and
PartitionedStateFluxes for the surface flux calculation.

#### References: Eisenman sea ice formulation
[Semtner 1976](https://journals.ametsoc.org/view/journals/phoc/6/3/1520-0485_1976_006_0379_amfttg_2_0_co_2.xml): A Model for the Thermodynamic Growth of Sea Ice in Numerical Investigations of Climate

[Eisenman & Wettlaufer 2009](https://eisenman.ucsd.edu/papers/Eisenman-Wettlaufer-2009.pdf): Nonlinear threshold behavior during the loss of Arctic sea ice

[Zhang et al. 2021](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JC016686): Sea Ice Properties in High-Resoluation Sea Ice Models
