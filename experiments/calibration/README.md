# ClimaCoupler Calibration Experiments

The `amip` experiment calibrates parameters of the coupled AMIP model against satellite
observations. A config file builds one `CalibrationTools.CalibrateConfig` named
`CALIBRATE_CONFIG`, which holds the whole experiment: the model yml, the loss fields, the
dates, the priors and the run settings. The loss fields are one line each:

```julia
observables = [
    (name = "swcre", reduce = seasonal_average ∘ zonal_average, model_error_fraction = 0.3),
    (name = "rtmt", reduce = seasonal_average ∘ global_mean, noise_sd = 0.2),
],
```

`reduce` is any function from an `OutputVar` to an `OutputVar`, applied to the observations and
the model output alike (`loss_values` in `preprocessing.jl`): `zonal_average`, `global_mean`,
your own, or a composition such as `seasonal_average ∘ zonal_average` for seasonal means. Maps take `model_error_fraction`, the model error
as a fraction of the observed mean (`noise_model.jl`). A global mean takes `noise_sd`, which sets
its weight against the maps. Comment out a line to drop a field. Observations and model output
stay in physical units. See the `CalibrateConfig` docstring for every setting.

To calibrate on a new field, add its observation's loader to `default_data_loader` and, if the
model does not output it directly in the same units, a `model_var` method, both in
`src/CalibrationTools.jl`.

For how to choose observables, noise and parameters, and what goes wrong, see
`amip/CALIBRATION_GUIDE.md`.

## Configs:

- `son2010_cre_toa.jl`: the default. swcre and lwcre zonal means plus the global net
  TOA flux (`rtmt`), SON 2010 seasonal means, 2 iterations, members from 2010-07-01 to
  2010-12-01, noise estimated over SON 2001-2010.
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
