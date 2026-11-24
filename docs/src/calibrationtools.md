```@meta
CurrentModule = CalibrationTools
```

# CalibrationTools

This module contains utilities to help with setting up calibration experiments
of the coupled model. These utilities are meant to be broadly useful and are not
specific to one particular calibration experiment.

!!! note "Other helpful resources"
    For calibration, other resources that you may find helpful are the
    documentation for
    [EnsembleKalmanProcesses](https://clima.github.io/EnsembleKalmanProcesses.jl/stable/),
    [ClimaCalibrate](https://clima.github.io/ClimaCalibrate.jl/stable/),
    [ClimaAnalysis](https://clima.github.io/ClimaAnalysis.jl/stable/), and the
    [ClimaCoupler calibration experiments](https://github.com/CliMA/ClimaCoupler.jl/tree/main/experiments/calibration).

## Data loading

!!! note "OutputVar"
    This section assumes you are familiar with `ClimaAnalysis.OutputVar`.

`CalibrationTools` provides a variety of data loaders from different artifacts
in `ClimaArtifacts` for calibration. For example, the [`ERA5DataLoader`](@ref)
is a data loader loading preprocessed ERA5 data. This data loader automatically
applies preprocessing to make it convenient to use for calibration. See the
documentation for [`ERA5DataLoader`](@ref) for details on the preprocessing
steps applied.

You can retrieve a variable with [`get`](@ref) and get a set of all available
preprocessed variables with [`available_vars`](@ref).

```@repl data_loader_example
import ClimaCoupler: CalibrationTools
data_loader = CalibrationTools.ERA5DataLoader()
CalibrationTools.available_vars(data_loader)
var = get(data_loader, "rsus");
```

In the example, we retrieve a `OutputVar` with the short name `rsus` which
represents the mean surface upward short-wave radiation flux.

!!! note "Other data loaders"
    If you want a data loader for other data sources, then please open an issue
    for it!

### I want to add a new variable to an existing data loader

To add a new variable, you must
1. Define a mapping between the ERA5 name and CliMA name,
2. Define a `preprocess` function for this variable.

To determine which variables are already available, refer to the artifact's
documentation. For [`ERA5DataLoader`](@ref) we can load the variable
representing mean evaporation rate or `mer` from the data source. We also want
to give it the name `er`. For step 1, we add `"mer" => "er"` as a mapping for
the data loader to recognize.

```@example add_variable
import ClimaCoupler: CalibrationTools
data_loader = CalibrationTools.ERA5DataLoader()
# ERA5_TO_CLIMA_NAMES define the existing pairings for the data loader
era5_to_clima_names = [CalibrationTools.ERA5_TO_CLIMA_NAMES..., "mer" => "er"]
data_loader = CalibrationTools.ERA5DataLoader(; era5_to_clima_names)
```

For the second step, we define a preprocessing function specific to the
variable.

!!! note "Preprocessing functions"
    See [ClimaAnalysis documentation](https://clima.github.io/ClimaAnalysis.jl/stable/)
    for available transformations on `OutputVar`s.

In our example, no preprocessing is applied.

```@example add_variable
CalibrationTools.preprocess(::CalibrationTools.ERA5DataLoader, var, ::Val{:er}) = var
nothing # hide
```

Now, you can use [`get`](@ref) to retrieve the `OutputVar` with the short name
"mer"`.

```@example add_variable
data_loader = CalibrationTools.ERA5DataLoader(; era5_to_clima_names)
get(data_loader, "er")
nothing # hide
```

## Building a loss

A calibration's loss is the list of `observables` in its [`CalibrateConfig`](@ref). Each
[`Observable`](@ref) names a field, a `reduce` function applied identically to the
observation and to the model output, and its error scale. The example
`experiments/calibration/amip/config/year2010_rsut_rlut_pr.jl` uses everything below.

### Reductions

The functions are in `experiments/calibration/amip/preprocessing.jl`.

| `reduce`             | Each value is                                                  |
|----------------------|----------------------------------------------------------------|
| `identity`           | a monthly mean, as the model writes them                       |
| `seasonal_average`   | a seasonal mean (DJF, MAM, JJA, SON), dated by its first month |
| `annual_average`     | a December-to-November mean, dated by its December             |
| `global_mean`        | the area-weighted global mean                                  |
| `zonal_average`      | the zonal mean                                                 |
| `remove_global_mean` | a map minus its own global mean                                |
| `annual_global_mean` | the global annual mean of a December-to-November year          |
| `annual_pattern`     | a December-to-November annual map minus its global mean        |
| `seasonal_anomaly`   | each season's map minus its year's mean                        |

`annual_global_mean`, `annual_pattern`, and `seasonal_anomaly` split a year of seasonal maps
into three orthogonal parts that add up to the seasonal maps and do not overlap, so a
calibration on all three cannot trade a global mean bias against a spatial pattern or
seasonal cycle. The helper `mean_pattern_seasonal_observables(field; global_mean_sd,
annual_pattern_fraction, seasonal_cycle_fraction)` builds all three `Observable`s for
`field`:

```julia
observables = [
    mean_pattern_seasonal_observables(
        "rsut";
        global_mean_sd = 12.6,
        annual_pattern_fraction = 0.42,
        seasonal_cycle_fraction = 0.10,
    )
    mean_pattern_seasonal_observables(
        "rlut";
        global_mean_sd = 4.8,
        annual_pattern_fraction = 0.187,
        seasonal_cycle_fraction = 0.187,
    )
    mean_pattern_seasonal_observables(
        "pr";
        global_mean_sd = 0.85,
        annual_pattern_fraction = 0.316,
        seasonal_cycle_fraction = 0.534,
    )
]
```

### Error scales and block weighting

Each observable has its own noise covariance block, estimated from the observation's
interannual spread over `covariance_date_ranges` plus a structural model-error term:

- `model_error_fraction` estimates an `SVDplusD` covariance across `covariance_date_ranges`
  with a structural model error equal to that fraction of the reduced observable's time
  mean. For a map whose sample mean crosses zero (such as `annual_pattern` or
  `seasonal_anomaly`), `generate_observations.jl` automatically floors the squared mean at
  its spatial median so zero crossings do not dominate the loss.
- `noise_sd` uses a diagonal covariance with per-entry standard deviation `noise_sd` in
  physical units.
- Setting `effective_entries` on an [`Observable`](@ref) multiplies that observable's noise
  covariance by `entries / effective_entries` when `generate_observations.jl` builds the
  observation vector (via `SampleBuilder.PerCollectionWeighting` for `SVDplusD` or directly
  on `ScalarDiagonal` for `DiagonalCovariance`). Dividing by `entries` normalizes any
  observable (whether a 1-entry global mean, a 1-map annual pattern, or a 4-map seasonal
  cycle) to its mean squared whitened residual per entry, and multiplying by
  `effective_entries` weights that normalized block as `effective_entries` independent
  unit-variance observations in both the loss and `EffectiveEntriesController` in
  `experiments/calibration/amip/step_size_rule.jl`.
- `mean_pattern_seasonal_observables` defaults to `effective_entries = 29` on all three
  parts (giving the global mean, annual pattern, and seasonal cycle equal 1/3 block weight
  and counting $29 + 29 + 29 = 87$ per field in `EffectiveEntriesController`, so
  `global_mean_sd` is given directly as the per-entry baseline error scale without manual
  `1 / sqrt(29)` rescaling).

## CalibrationTools API

```@docs
CalibrationTools.CalibrateConfig
CalibrationTools.CalibrateConfig()
CalibrationTools.Observable
CalibrationTools.Observable()
CalibrationTools.ERA5DataLoader
CalibrationTools.ERA5DataLoader()
CalibrationTools.CERESDataLoader
CalibrationTools.CERESDataLoader()
CalibrationTools.GPCPDataLoader
CalibrationTools.GPCPDataLoader()
CalibrationTools.ERA5PressureLevelDataLoader
CalibrationTools.ERA5PressureLevelDataLoader()
CalibrationTools.ModisDataLoader
CalibrationTools.ModisDataLoader()
CalibrationTools.CompositeDataLoader
CalibrationTools.CompositeDataLoader()
CalibrationTools.find_source_loader
CalibrationTools.available_vars
CalibrationTools.get
CalibrationTools.update_timespan!
CalibrationTools.add_parameter_filepath!
```
