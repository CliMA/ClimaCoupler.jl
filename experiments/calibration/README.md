# ClimaCoupler Calibration Experiments

The `amip` experiment calibrates parameters of the coupled AMIP model against
satellite observations. A config file builds one
`CalibrationTools.CalibrateConfig` named `CALIBRATE_CONFIG`, which holds the
whole experiment: the model YAML file, the loss fields, the dates, the priors and the run settings. The loss fields are one `CalibrationTools.Observable` each:

```julia
observables = [
    Observable("swcre", seasonal_average ∘ zonal_average; model_error_fraction = 0.3),
    Observable("rtmt", seasonal_average ∘ global_mean; noise_sd = 0.2),
]
```

`reduce` is any function from an `OutputVar` to an `OutputVar`, applied to the
observations and the model output alike (`mask_and_reduce` in `preprocessing.jl`):
`zonal_average`, `global_mean`, your own, or a composition such as
`seasonal_average ∘ zonal_average` for seasonal means. `model_error_fraction` estimates an
`SVDplusD` covariance across `covariance_date_ranges` with model error as a fraction of the
reduced observable's mean (`make_observation_vector` in `generate_observations.jl`).
`noise_sd` gives diagonal noise with standard deviation `noise_sd` in physical units.
Comment out an entry to drop a field. Observations and model output stay in physical units.
See the `Observable` and `CalibrateConfig` docstrings for every setting.

To calibrate on a new field, add its observation's loader to
`default_data_loader` and, if the model does not output it directly in the same
units, a `model_var` method, both in `src/CalibrationTools.jl`.

## Configs:

- `year2010_rsut_rlut_pr.jl`: the default. rsut, rlut, and pr from December 2009 to
  November 2010, each as its global annual mean, its annual map minus that mean,
  and each season's map minus the annual map, with each of the nine parts weighted
  as 29 effective entries. Members from 2009-11-01 to 2010-12-01, noise
  estimated over the 23 December-to-November years from December 2000 to
  November 2023.
- `pipeline_test.jl`: an end-to-end test that runs a single iteration using `ta`
  and `hur` observations of October 2010. To use this config, you
  should set the environment variable `TEST_CALIBRATION` to anything before
  generating the observations and starting the calibration. It runs the
  ensemble members one after another on the CPU in the current process, so it
  needs no workers or job submissions.

Set `CALIBRATION_CONFIG` in the environment to run another config file without editing
`run_calibration.jl` or `generate_observations.jl`.

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
