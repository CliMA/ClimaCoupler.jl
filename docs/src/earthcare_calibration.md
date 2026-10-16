# EarthCARE cloud calibration

A staged calibration of the 1M prognostic-EDMF AMIP configuration against EarthCARE
radar-lidar cloud retrievals. The stages add observables and parameters one step at a
time, so each stage can be checked before the next is run.

## Targets

EarthCARE ACM_CLP (JAXA CPR-ATLID-MSI synergy) granules from December 2025 to May 2026
are re-binned from raw files onto the model grid with fill values treated as zero. The
re-bin script and products live outside the repository for now
(`~/earthcare_eval` on the clima cluster):

- `targets/earthcare_targets.nc`: frozen water path `fwp` (kg m^-2, the column integral
  of radar-lidar ice water content, which includes falling snow) and cloud fraction at
  1, 3 and 8 km (`cl1km`, `cl3km`, `cl8km`, fraction 0 to 1, from the CPR-ATLID-MSI mask).
  The time axis is relabelled to 2010 so each month pairs with a 2010 model run.
- `targets/era5_fwp_2001_2010.nc`: ERA5 cloud ice plus snow water path for the same
  calendar months over 2001 to 2010, used only to estimate the noise covariance.
  ERA5 and EarthCARE agree on this quantity to within about 10 percent globally.

The existing gridded product in the EarthCARE archive (`binned/`) is not used: it keeps
only pixels with both phases retrieved and does not mask the fill sentinel.

## Matching model quantities

- `fwp` is `iwp + swp`: the radar sees all frozen water in the volume, so the model's
  cloud ice must be added to its snow. This needs the 1M microphysics.
- Cloud fraction is `cl` interpolated to height. The EarthCARE mask is a volume
  occurrence per 100 m bin, so expect it to sit above the model's grid-mean fraction.
- EarthCARE liquid water is only retrieved where the lidar is unattenuated, so neither
  its liquid water content nor its volume liquid fraction is a target. A cloud-top
  liquid fraction will need a cloud-top diagnostic in ClimaAtmos (stage 3).

## Stages

1. `config/earthcare_fwp.jl`: `fwp`, `swcre`, `lwcre`, `pr` against
   `snow_autoconversion_timescale`, `cloud_ice_snow_collision_efficiency`,
   `snow_apparent_density`, `rain_autoconversion_timescale`. Snow falls with the Chen2022
   speed by default on ClimaAtmos main, and its apparent density sets those coefficients;
   the fixed fall speed and the aspect ratio are not read in this configuration (a first
   two-iteration run with `fixed_snow_terminal_velocity` showed bit-identical members).
   Targets alternate April and January 2010, nine ensemble members.
2. Add cloud fraction at three heights and MAC `lwp`; add `cloud_fraction_eps_rel`,
   `entr_coeff`, `detr_massflux_vertdiv_coeff`, `cloud_liquid_rain_collision_efficiency`.
3. Add cloud-top liquid fraction versus temperature and ERA5 `ta`, `hur`; add ice
   nucleation and fall-speed parameters.

## Running stage 1

```bash
export CALIBRATION_CONFIG=earthcare_fwp.jl
julia --project=experiments/AMIP experiments/calibration/amip/generate_observations.jl
julia --project=experiments/AMIP experiments/calibration/amip/run_calibration.jl
```

`EARTHCARE_TARGETS_NC` and `ERA5_FWP_NC` override the target and covariance file paths;
`CALIBRATION_OUTPUT_DIR` overrides the output directory, which defaults to local scratch.

## Noise model

`generate_observations.jl` accepts `COVARIANCE_DATA_LOADERS` to estimate the covariance
of a variable from a different source than its target, and `COVARIANCE_BY_MONTH` to
estimate one covariance per calendar month. Stage 1 uses both: the EarthCARE record has
one sample per month, so the interannual spread comes from ERA5.

## Initial conditions

The members start seven days before each target month, so on 25 December 2009 and
25 March 2010. `initial_condition: "WeatherModel"` reads
`era5_init_processed_internal_<date>_0000.nc` from the `wxquest_initial_conditions`
artifact, which does not hold those dates. `config/amip_configs/amip_calibration.yml`
therefore sets `era5_initial_condition_dir` to a directory on local scratch (the home
directory is NFS with a quota, and a 9 GB write there failed); a missing file is
fetched from the Copernicus Climate Data Store by ClimaInitialConditions, which needs a
`~/.cdsapirc`. Fetch both dates once before launching, so the members do not queue on
the download: `ClimaInitialConditions.ERA5.fetch_initial_conditions(date; dir)`.

The driver is a long CPU batch job (`~/earthcare_eval/run_stage1.sbatch` on the clima
cluster); ClimaCalibrate submits each member as its own GPU job.
