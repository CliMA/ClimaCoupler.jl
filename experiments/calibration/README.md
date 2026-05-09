# ClimaCoupler Calibration Experiments

The `amip` experiment calibrates parameters against ERA5 pressure-level
observations of `ta` and `hur` at 200, 500, and 850 hPa with a latitude-weighted
scalar covariance matrix. The observations are z-score normalized for each variable and pressure level.

## Configs:

- `pressure_levels.jl`: 6 iterations using `ta` and `hur` observations of
  October 2010.
- `pipeline_test.jl`: an end-to-end test that runs a single iteration using `ta`
  and `hur` observations of October 2010. To use this config, you
  should set the environment variable `TEST_CALIBRATION` to anything before
  generating the observations and starting the calibration.

## Running a calibration

### Generate observations

Before running a calibration, you need to generate the observations with:


```bash
julia --project=experiments/AMIP experiments/calibration/amip/generate_observations.jl
```

### Run calibration

Then, you can run the calibration with

```bash
julia --project=experiments/AMIP experiments/calibration/amip/run_calibration.jl
```
