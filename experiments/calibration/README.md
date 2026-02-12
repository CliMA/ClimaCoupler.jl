# ClimaCoupler Calibration Experiments


This folder contains pipelines to reproduce coupled calibration experiments.
Each pipeline has its own subfolder:

- perfect_model: A trivial perfect-model calibration of the atmosphere coupled
with the bucket model. The calibration uses 30-day and lat/lon averages of
top-of-atmosphere shortwave radiation to calibrate the `total_solar_irradiance`
parameter in a perfect model setting. The current run script uses the
 `ClimaCalibrate.SlurmManager` to add Slurm workers which run each ensemble
 member in parallel.
- subseasonal: Calibrates the inverse entrainment timescale to ERA5 October monthly surface fluxes and surface temperature from 2018 to 2024.

To run a pipeline on a Slurm cluster, ensure that the given runscript
`experiments/calibration/<pipeline>/run_calibration.jl` is configured for your
cluster and run:
```
julia --project=experiments/ClimaEarth experiments/calibration/<pipeline>/run_calibration.jl
```
