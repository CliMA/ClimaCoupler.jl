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
    [ClimaCalibrate](https://clima.github.io/ClimaCalibrate.jl/dev/),
    [ClimaAnalysis](https://clima.github.io/ClimaAnalysis.jl/stable/), and the
    [ClimaCoupler calibration experiments](https://github.com/CliMA/ClimaCoupler.jl/tree/main/experiments/calibration).

## Data loading

!!! note "OutputVar"
    This section assumes you are familiar with `ClimaAnalysis.OutputVar`.

As of now, `CalibrationTools` provide a single data loader which is the
[`ERA5DataLoader`](@ref) for loading preprocessed ERA5 data. This data loader
automatically applies preprocessing to make it convenient to use for
calibration. See the documentation for [`ERA5DataLoader`](@ref) for details on
the preprocessing steps applied.

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
    See [ClimaAnalysis documentation](https://clima.github.io/ClimaAnalysis.jl/dev/)
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

## CalibrationTools API

```@docs
CalibrationTools.CalibrateConfig
CalibrationTools.CalibrateConfig()
CalibrationTools.ERA5DataLoader
CalibrationTools.ERA5DataLoader()
CalibrationTools.available_vars
CalibrationTools.get(loader::CalibrationTools.ERA5DataLoader, short_name::String)
```
