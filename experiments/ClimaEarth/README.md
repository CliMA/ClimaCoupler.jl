# ClimaEarth Experiments

The main script in the ClimaEarth directory is [run_amip.jl](run_amip.jl). This script can be used to
run the AMIP configuration, a simplified slabplanet configuration, or a CMIP configuration
using Oceananigans.jl. All of these options use sequential coupling.

The configuration to run is determined by the input `config_file` supplied to the
run_amip.jl driver, specifically by the option `mode_name`. Please see the top-level
[README.md](https://github.com/CliMA/ClimaCoupler.jl?tab=readme-ov-file#climacouplerjl)
for information about how to run the driver.

Each configuration is explained here.

## CMIP components
CMIP is currently the most complex configuration of the ClimaEarth model.
It runs a [ClimaAtmos.jl](https://github.com/CliMA/ClimaAtmos.jl) atmosphere model,
[ClimaLand.jl](https://github.com/CliMA/ClimaLand.jl) bucket or integrated land model,
[Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) ocean model,
and either [ClimaSeaIce.jl](https://github.com/CliMA/ClimaSeaIce.jl) or a simple thermal sea ice model.
Please find more details about each model below.

!!! note Behavior at the poles with an ocean model on a capped latitude/longitude grid
    The Oceananigans and ClimaSeaIce models currently run on a capped
    latitude-longitude grid, which spans from 80°S to 80°N. To avoid having a gap
    in surface models at the poles, we fill the poles with the selected land model
    when running with the Oceananigans model. As a result, the land model cannot be
    started from saved initial conditions when run in this configuration.
    This will change in the near future when we switch to use a tripolar grid
    for the Oceananigans and ClimaSeaIce models.

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

#### References: Atmosphere formulation

[Manabe 1969](https://journals.ametsoc.org/view/journals/mwre/97/11/1520-0493_1969_097_0739_catoc_2_3_co_2.xml): Climate and ocean circulation. I. The atmospheric circulation and the hydrology of the earth's surface.

[Pincus et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621): Balancing Accuracy, Efficiency, and Flexibility in Radiation Calculations for Dynamical Models

[Tan et al. 2018](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017MS001162): An Extended Eddy-Diffusivity Mass-Flux Scheme for Unified Representation of Subgrid-Scale Turbulence and Convection

[Cohen et al. 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002162): Unified Entrainment and Detrainment Closures for Extended Eddy-Diffusivity Mass-Flux Schemes

[Lopez-Gomez et al. 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002161): A Generalized Mixing Length Closure for Eddy-Diffusivity Mass-Flux Schemes of Turbulence and Convection

### Land

#### Bucket option
Dynamical core:
- Equation: Following Manabe bucket hydrology scheme (Manabe 1969)
- Prognostic variables: Temperature, water content, snow water content
- Spatial discretization: Finite difference, multiple independent columns on the sphere
- Time stepping: Explicit additive Runge–Kutta method

Land surface albedo:
- Bare ground: Prescribed from files
- Snow: Constant

#### Integrated land option
This is a more complex land model than the bucket, with multiple models and
parameterizations nested in a modular structure. At the top level, the land model
consists of soil, canopy, and snow models.

Dynamical core:
- Equation: See the ClimaLand.jl documentation section ["Model Equations"](https://clima.github.io/ClimaLand.jl/stable/)
- Prognostic variables:
  - Soil: internal energy, water content, ice content, carbon content
  - Canopy: temperature, water content,
  - Snow: snow water equivalent, snow liquid water, energy
- Spatial discretization: Finite difference, multiple independent columns on the sphere
- Time stepping: Mixed implicit/explicit (IMEX) additive Runge–Kutta method

Please see the [ClimaLand.jl documentation](https://clima.github.io/ClimaLand.jl/stable/)
for more information and examples using combinations of ClimaLand
models/parameterizations.

### Ocean
Oceananigans.jl is used for the CMIP ocean model. Please see that package's
[documentation](https://clima.github.io/OceananigansDocumentation/stable/)
for details about the physics of the model.

### Sea Ice

#### ClimaSeaIce
ClimaSeaIce.jl is the more complex of the two sea ice options used in the CMIP configuration.

#### Prescribed sea ice
Thermodynamics: 0-layer model, with prognostic ice surface temperature and fixed ice thickness

Sea ice concentration is prescribed by the `SIC` input file and changes with time.
This is provided as a simpler alternative to ClimaSeaIce.

## AMIP components
AMIP runs a ClimaAtmos.jl atmosphere model, ClimaLand.jl bucket land model,
a prescribed ocean model, and a simple thermal sea ice model.

Note that this is very similar to the CMIP setup, with the replacement of a simpler
ocean model, and the restriction that the prescribed sea ice model must be used.

### Ocean (prescribed SST)
The ocean model has no prognostic variables. It stores sea surface temperature, which is
read from the `SST` input file and changes over time.

## Slabplanet components
The Slabplanet configuration is more idealized than the AMIP configuration, but provides
valuable insight about conservation and individual model behavior.

### Atmosphere
The atmosphere model used in the slabplanet configuration is the same as the one used for CMIP;
please see that section for details.

### Land
The land model used in the slabplanet configuration is the bucket described in the CMIP section;
please see that section for details. Note that we could also use the integrated land without requiring
changes to the code, but this is not currently being exercised.

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

## Subseasonal mode

Generates 3-4 week forecasts from ERA5-derived initial conditions. Similar to `amip` setup. Given a directory containing initial condition files, the filenames associated with current date are inferred from `start_date`. Initial condition files are generated by `https://github.com/CliMA/WeatherQuest`.

- `era5_initial_condition_dir`: directory that contains files named like:
  - `sst_processed_YYYYMMDD_0000.nc`
  - `sic_processed_YYYYMMDD_0000.nc`
  - `era5_land_processed_YYYYMMDD_0000.nc`
  - `era5_bucket_processed_YYYYMMDD_0000.nc`

- **How filename inference works**: given `start_date` formatted as `YYYYMMDD`, the model constructs paths
  - `joinpath(era5_initial_condition_dir, "sst_processed_YYYYMMDD_0000.nc")` for the initial conditions.

- **Contents and expected variables**:
  - `sst_processed_YYYYMMDD_0000.nc`:
    - Variable: `SST`
    - Units: degrees C in file; converted internally to K
  - `sic_processed_YYYYMMDD_0000.nc`:
    - Variable: `SEAICE`
    - Units: percent (0–100) in file; converted internally to fraction (0–1)
  - `era5_land_processed_YYYYMMDD_0000.nc` (for integrated land model):
    - Variables used directly: `skt` (skin temperature, K), `tsn` (snow surface temperature, K)
    - Required fields for land initialization:
      - `swe` (m): snow water equivalent
      - `swvl` (m^3/m^3): volumetric fraction of liquid water, dims `(z, lat, lon)`
      - `si` (m^3/m^3): volumetric fraction of ice, dims `(z, lat, lon)`
      - `sie` (J/m^3): soil volumetric internal energy, dims `(z, lat, lon)`
      - `stl` (K): soil temperature, dims `(z, lat, lon)`
      - `tsn` (K): temperature of snow layer
      - `skt` (K): skin temperature
    - era5 land level midpoints: 0.035, 0.175, 0.64, 1.945 (m)

  - `era5_bucket_processed_YYYYMMDD_0000.nc` (for bucket land model):
    - Required fields for initialization:
      - `W` (m): subsurface water content, dims `(lat, lon)`
      - `Ws` (m): surface water content, dims `(lat, lon)`
      - `S` (m): snow water equivalent, dims `(lat, lon)`
      - `T` (K): soil temperature profile, dims `(z, lat, lon)`
      - `tsn` (K): temperature of snow layer, dims `(lat, lon)`
      - `skt` (K): skin temperature, dims `(lat, lon)`
    - era5 land level midpoints: 0.035, 0.175, 0.64, 1.945 (m)

## Configuration files
We use configuration files to specify all simulation parameters and options. The configuration files are organized hierarchically:

- **Atmosphere configurations**: Located in `config/atmos_configs/`, these files contain atmosphere-specific settings
- **Coupler configurations**: Located in other `config/` subdirectories, these files reference atmosphere configurations and add settings for surface models and the coupling system.

The coupler configuration takes precedence over any conflicting settings from the atmosphere configuration.
