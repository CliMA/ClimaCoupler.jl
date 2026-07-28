# ERA5 initial conditions

Subseasonal (WeatherModel) simulations start from ERA5 reanalysis data. The
`Era5InitialConditions` module fetches this data on demand from the
[Copernicus Climate Data Store](https://cds.climate.copernicus.eu) (CDS) for
the run's `start_date`, preprocesses it into the files the coupler and its
component models read, validates it, and caches it locally. Repeat runs for
the same date make no network requests.

## Credentials

The CDS fetch needs a (free) CDS account:

1. Register at [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu)
   and accept the ERA5 licence
   ("Licence to use Copernicus Products").
2. Create `~/.cdsapirc` with your personal access token, following the
   [CDS API instructions](https://cds.climate.copernicus.eu/how-to-api):

   ```
   url: https://cds.climate.copernicus.eu/api
   key: <your-personal-access-token>
   ```

   Alternatively, set the `CDSAPI_URL` and `CDSAPI_KEY` environment
   variables.

CDS requests are queued on the Copernicus servers. A fetch typically takes
minutes, but can take hours when the queue is busy. The download is roughly
1-2 GB per date.

## Resolution order

For subseasonal runs, the initial condition directory is resolved in this
order:

1. An explicit `--era5_initial_condition_dir` (used as-is, no fetch).
2. The local cache, or a CDS fetch when the date is not cached.
3. The `wxquest_initial_conditions` ClimaArtifact (a fixed set of
   pre-generated dates).

## Caching

Fetched files land in a per-package
[Scratch.jl](https://github.com/JuliaPackaging/Scratch.jl) directory and
survive across runs. Set the `CLIMACOUPLER_ERA5_CACHE_DIR` environment
variable to use a different directory, for example a shared cache on a
cluster. Files for many dates coexist in one cache directory.

A fetch downloads into a temporary directory and moves the files into the
cache only after validation, so an interrupted run never leaves a partial
cache. Pass `force = true` to
[`Era5InitialConditions.fetch_era5_initial_conditions`](@ref) to delete and
refetch a date.

## Produced files

For a start date `YYYYMMDD`, the fetch produces:

- `era5_raw_YYYYMMDD_0000.nc`: pressure-level atmosphere state (`u`, `v`,
  `w`, `t`, `q`, `z`, cloud condensate) plus `skt`, `sp`, and
  `surface_geopotential`, read by ClimaAtmos's WeatherModel initial
  condition. ClimaAtmos interpolates it to its grid, including the vertical
  interpolation from pressure levels.
- `sst_processed_YYYYMMDD_0000.nc`: `SST` in Celsius, land filled by nearest
  neighbor.
- `sic_processed_YYYYMMDD_0000.nc`: `SEAICE` in percent.
- `era5_land_processed_YYYYMMDD_0000.nc`: integrated-land initial conditions
  (`skt`, `tsn`, `swe`, `swvl`, `stl`, `si`, `sie`, `lai`); ocean points
  are 0.
- `era5_bucket_processed_YYYYMMDD_0000.nc`: bucket initial conditions (`W`,
  `Ws`, `S`, `T`, `tsn`, `skt`).
- `albedo_processed_YYYYMMDD_0000.nc`: `sw_alb_clr` (ERA5 forecast albedo),
  used when `bucket_albedo_type: "era5"`.

The preprocessing follows the WeatherQuest pipeline that produced the
`wxquest_initial_conditions` artifact, with one difference: the atmosphere
state uses the 37 ERA5 pressure levels instead of the 137 model levels, and
ClimaAtmos performs the vertical interpolation.

## API

```@docs
ClimaCoupler.Era5InitialConditions.fetch_era5_initial_conditions
ClimaCoupler.Era5InitialConditions.default_cache_dir
ClimaCoupler.Era5InitialConditions.cds_credentials_available
ClimaCoupler.Era5InitialConditions.validate_era5_dir
```
