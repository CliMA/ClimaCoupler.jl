# Available Simulation Types

The simulation type or "mode" controls which component models are used and how they are configured.
It is selected via the `mode_name` configuration option (see [Input](@ref)).

| `mode_name` | Description |
|---|---|
| `amip` | Atmosphere + prescribed ocean/sea ice + land |
| `cmip` | Atmosphere + dynamic ocean/sea ice + land |
| `slabplanet` | Atmosphere + slab ocean + land |
| `slabplanet_aqua` | Atmosphere + slab ocean only (aquaplanet) |
| `slabplanet_terra` | Atmosphere + land only |
| `subseasonal` | Short-range forecast from ERA5 initial conditions using AMIP setup|

All simulation types use [ClimaAtmos.jl](https://github.com/CliMA/ClimaAtmos.jl) as the
atmosphere model. For AMIP, CMIP, and subseasonal runs the land model may be either the
bucket or integrated land model from [ClimaLand.jl](https://github.com/CliMA/ClimaLand.jl);
for slabplanet, aquaplanet, and terraplanet runs the bucket land model must be used.
See [Available component models](@ref) for details on each component model.

## CMIP (`mode_name: "cmip"`)

CMIP (Coupled Model Intercomparison Project) is the most complex simulation type
supported by ClimaCoupler.jl. In addition to the prognostic atmosphere and land
models, the ocean evolves prognostically in response to atmospheric forcing,
and sea ice is thermodynamically active.

**Component models:**
- Atmosphere: `ClimaAtmosSimulation`
- Land: `BucketSimulation` or `ClimaLandSimulation` (controlled by `land_model`)
- Ocean: `OceananigansSimulation`
- Sea ice: `ClimaSeaIceSimulation`

!!! tip "GPU recommended"
    The CMIP configuration is computationally expensive due to the complexity of all
    component models. Running on a GPU is strongly recommended; see
    the `device` option in [Input](@ref) for how to select the compute device.

!!! note "Behavior at the poles"
    The Oceananigans and ClimaSeaIce models currently run on a capped latitude-longitude
    grid spanning 80°S to 80°N. To avoid a gap at the poles, the selected land model is
    used to fill the polar regions. As a result, the land model cannot be started from saved
    initial conditions in this configuration. This will change in the future when the models
    switch to a tripolar grid.

## AMIP (`mode_name: "amip"`)

AMIP (Atmospheric Model Intercomparison Project) is a standard experimental protocol of the
[Program for Climate Model Diagnosis & Intercomparison (PCMDI)](https://pcmdi.llnl.gov/).
It is used to evaluate atmosphere and land model components while sea surface temperatures
(SST) and sea ice concentration (SIC) are prescribed from observational data (e.g., HadISST).

**Component models:**
- Atmosphere: `ClimaAtmosSimulation`
- Land: `BucketSimulation` or `ClimaLandSimulation` (controlled by `land_model`)
- Ocean: `PrescribedOceanSimulation`
- Sea ice: `PrescribedIceSimulation`

## Slabplanet (`mode_name: "slabplanet"`)

The slabplanet configuration is a more idealized setup than AMIP, designed for studying
conservation properties and individual model behavior. The ocean is a simple thermal slab
with a prognostic surface temperature but no dynamics, and there is no sea ice.

**Component models:**
- Atmosphere: `ClimaAtmosSimulation`
- Land: `BucketSimulation` or `ClimaLandSimulation`
- Ocean: `SlabOceanSimulation`
- Sea ice: none (ocean fills ice-covered regions)

### Slabplanet aqua (`mode_name: "slabplanet_aqua"`)

An aquaplanet variant: the slab ocean covers the entire surface with no land or sea ice.

**Component models:**
- Atmosphere: `ClimaAtmosSimulation`
- Ocean: `SlabOceanSimulation` (entire surface)
- Land: none
- Sea ice: none

### Slabplanet terra (`mode_name: "slabplanet_terra"`)

A land-only variant: the land model covers the entire surface with no ocean or sea ice.

**Component models:**
- Atmosphere: `ClimaAtmosSimulation`
- Land: `BucketSimulation` (entire surface)
- Ocean: none
- Sea ice: none

## Subseasonal (`mode_name: "subseasonal"`)

Generates 3–4 week forecasts initialized from ERA5 reanalysis data. The setup is otherwise
similar to AMIP, but uses specific ERA5-derived initial conditions for the land model. The
`era5_initial_condition_dir` option must point to a directory containing the initial
condition files described below.

Initial condition files can be generated using the
[WeatherQuest](https://github.com/CliMA/WeatherQuest) package.

### Expected input files

Given a `start_date` formatted as `YYYYMMDD`, the following files are expected in
`era5_initial_condition_dir`:

| File | Contents |
|---|---|
| `sst_processed_YYYYMMDD_0000.nc` | SST variable `SST` in °C; converted internally to K |
| `sic_processed_YYYYMMDD_0000.nc` | Sea ice concentration `SEAICE` in percent; converted to fraction |
| `era5_bucket_processed_YYYYMMDD_0000.nc` | Bucket land IC; auto-inferred if `bucket_initial_condition` not set |
| `era5_land_processed_YYYYMMDD_0000.nc` | Integrated land IC; required for `land_model: integrated` |
| `albedo_processed_YYYYMMDD_0000.nc` | Optional surface albedo; used when `bucket_albedo_type: era5` |

#### Bucket land IC (`era5_bucket_processed_YYYYMMDD_0000.nc`)

| Variable | Units | Dimensions | Description |
|---|---|---|---|
| `W` | m | `(lat, lon)` | Subsurface water content |
| `Ws` | m | `(lat, lon)` | Surface water content |
| `S` | m | `(lat, lon)` | Snow water equivalent |
| `T` | K | `(z, lat, lon)` | Soil temperature profile |
| `tsn` | K | `(lat, lon)` | Snow layer temperature |
| `skt` | K | `(lat, lon)` | Skin temperature |

ERA5 land level midpoints: 0.035, 0.175, 0.64, 1.945 m.

#### Integrated land IC (`era5_land_processed_YYYYMMDD_0000.nc`)

| Variable | Units | Dimensions | Description |
|---|---|---|---|
| `swe` | m | `(lat, lon)` | Snow water equivalent |
| `swvl` | m³/m³ | `(z, lat, lon)` | Volumetric liquid water fraction |
| `si` | m³/m³ | `(z, lat, lon)` | Volumetric ice fraction |
| `sie` | J/m³ | `(z, lat, lon)` | Soil volumetric internal energy |
| `stl` | K | `(z, lat, lon)` | Soil temperature |
| `tsn` | K | `(lat, lon)` | Snow layer temperature |
| `skt` | K | `(lat, lon)` | Skin temperature |

ERA5 land level midpoints: 0.035, 0.175, 0.64, 1.945 m.

#### Surface albedo (`albedo_processed_YYYYMMDD_0000.nc`)

Used when `bucket_albedo_type: era5`. Contains `sw_alb_clr` (clear-sky surface albedo,
fraction), with dimensions `(time, lat, lon)` representing monthly data that is
temporally interpolated.
