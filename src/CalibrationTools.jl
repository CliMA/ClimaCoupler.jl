"""
    CalibrationTools

This module contains utilities to simplify calibration experiments of a coupled
simulation.
"""
module CalibrationTools

import Dates

import ClimaAnalysis
import ClimaAnalysis: NCCatalog
import ClimaCoupler
import ClimaCoupler.Utilities
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaDiagnostics
import ClimaUtilities.TimeManager: ITime

# TODO: Remove this once the name is added to ClimaAnalysis
push!(ClimaAnalysis.Var.ALTITUDE_NAMES, "height")

"""
    struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period}

A configuration struct for keeping track of multiple fields that are of interest
to a user running calibration, or that are needed in multiple places (e.g., for
ensemble members and generating observations).
"""
struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period, OBS, PRIORS}
    "Configuration file to use for ClimaCoupler simulation"
    config_file::String

    "The loss fields, one named tuple each (see the constructor)"
    observables::OBS

    "Pressure levels (Pa) kept from fields on pressure levels"
    pressure_levels::Vector{Float64}

    "The number of iterations to run the calibration for"
    n_iterations::Int64

    "The target date range of each iteration (`n_iterations + 1` entries), which also sets
    the start and end dates of its simulations"
    sample_date_ranges::Vector{NTuple{2, Dates.DateTime}}

    "The date ranges whose observations estimate the noise covariance"
    covariance_date_ranges::Vector{NTuple{2, Dates.DateTime}}

    "The amount of time to run a simulation after the end of its target date range"
    extend::EXTEND

    "The amount of time to run a simulation before the start of its target date range"
    spinup::SPINUP

    "The directory to store the iterations and members of the calibration."
    output_dir::String

    "An integer value for ensuring calibrations are the same between multiple
    calibrations with the same settings"
    rng_seed::Int64

    "The prior distribution of the calibrated parameters"
    priors::PRIORS

    "Walltime of each worker job, in minutes"
    worker_walltime::Int64
end

"""
    check_priors(priors)

Error if any prior is a unit normal over its bounds. EKP's `constrained_gaussian` returns
its unit-normal starting point unchanged for parameters below about 1e-2, so the prior
silently ignores the requested mean and spread, and that exact unit normal is its
signature. `priors` is an EKP `ParameterDistribution`. It is inspected by its fields, so
ClimaCoupler does not depend on EKP.
"""
function check_priors(priors)
    for (name, d) in zip(priors.name, priors.distribution)
        hasproperty(d, :distribution) || continue
        dist = d.distribution
        nameof(typeof(dist)) == :Normal &&
            dist.μ == 0 &&
            dist.σ == 1 &&
            error("The prior of $name is a unit normal over its bounds, the sign of \
                  constrained_gaussian failing on a small parameter. Build it from the \
                  unconstrained normal directly.")
    end
    return nothing
end

"""
    CalibrateConfig(;
        config_file,
        observables = [],
        pressure_levels = 100.0 .* [200.0, 500.0, 850.0],
        n_iterations::Integer,
        sample_date_ranges,
        covariance_date_ranges = [],
        extend::Dates.Period,
        spinup::Dates.Period,
        output_dir,
        rng_seed = 42,
        priors = nothing,
        worker_walltime = 350,
    )

Initializes a `CalibrateConfig`, which holds everything a calibration experiment needs.

Keyword arguments
=====================

- `config_file`: Configuration file to use for ClimaCoupler simulation.

- `observables`: The loss fields, one named tuple each, for example

      (name = "swcre", reduce = zonal_average, model_error_fraction = 0.3)
      (name = "rtmt", reduce = global_mean, noise_sd = 0.2)

  `name` is the short name of the field. `reduce` is a function from an `OutputVar` to an
  `OutputVar`, applied identically to the observations and the model output. Compose it with
  `seasonal_average` to compare seasonal means, e.g. `seasonal_average ∘ zonal_average`,
  and set `extend` to match. A field that
  keeps latitude takes `model_error_fraction`, the model error as a fraction of the observed
  mean. A field reduced to one value per date takes `noise_sd`, independent noise in its own
  units, whose value sets its weight. Required to generate observations.

  A field's observation comes from `default_data_loader()` and its model output from
  `model_var`.

- `pressure_levels`: Pressure levels (Pa) kept from fields on pressure levels.

- `n_iterations`: The number of iterations to run the calibration for.

- `sample_date_ranges`: The target date ranges, each once. Iterations cycle through them in
  order, one per iteration. The dates should be the same as found in the time series data
  of the observations after the `reduce` functions, e.g. `(DateTime(2010, 9, 1),
  DateTime(2010, 9, 1))` for SON 2010 seasonal means.

- `covariance_date_ranges`: The date ranges whose observations estimate the noise
  covariance, typically the same season in about ten years. Required to generate
  observations.

- `extend`: The amount of time to run the simulation after the end date
  determined by `sample_date_ranges`. For seasonal averages, `extend` should be
  `Dates.Month(3)` and for monthly averages, `extend` should be
  `Dates.Month(1)`.

- `spinup`: The amount of time to run the simulation before the start date
  determined by `sample_date_ranges`.

- `output_dir`: The location to save the calibration at.

- `rng_seed`: An integer to ensure that calibration runs with the same settings
  are the same.

- `priors`: The prior distribution of the calibrated parameters, as built by
  `EnsembleKalmanProcesses.combine_distributions`. Required to run a calibration. Each
  prior is checked for the unit normal `constrained_gaussian` returns when it fails on a
  small parameter.

- `worker_walltime`: Walltime of each worker job, in minutes. It must cover
  `n_iterations` member runs.
"""
function CalibrateConfig(;
    config_file,
    observables = [],
    pressure_levels = 100.0 .* [200.0, 500.0, 850.0],
    n_iterations::Integer,
    sample_date_ranges,
    covariance_date_ranges = [],
    extend::Dates.Period,
    spinup::Dates.Period,
    output_dir,
    rng_seed = 42,
    priors = nothing,
    worker_walltime = 350,
)
    isfile(config_file) || error("Configuration file ($config_file) does not exist")
    isempty(sample_date_ranges) &&
        error("Cannot run calibration with no date ranges for the samples")

    sign(extend) == -1 && error("The period to extend ($extend) by should be positive")
    sign(spinup) == -1 && error("The period to spin up ($spinup) by should be positive")

    n_iterations > 0 || error("The number of iterations ($n_iterations) should be positive")

    # Typed, so that an empty list (no covariance dates) still matches the field type.
    to_datetimes(ranges) = NTuple{2, Dates.DateTime}[
        (Dates.DateTime(first(r)), Dates.DateTime(last(r))) for r in ranges
    ]
    targets = to_datetimes(sample_date_ranges)
    # The calibration reads one entry per iteration plus one.
    sample_date_ranges = [targets[mod1(i, length(targets))] for i in 1:(n_iterations + 1)]
    covariance_date_ranges = to_datetimes(covariance_date_ranges)

    for (start_date, stop_date) in vcat(sample_date_ranges, covariance_date_ranges)
        start_date <= stop_date || error(
            "The start date ($start_date) should be before the stop date ($stop_date)",
        )
    end

    for o in observables
        haskey(o, :name) && haskey(o, :reduce) ||
            error("Observable $o needs a name and a reduce function")
        haskey(o, :model_error_fraction) == haskey(o, :noise_sd) && error(
            "Observable $(o.name) needs exactly one of model_error_fraction and noise_sd",
        )
    end
    names = [o.name for o in observables]
    allunique(names) || error("Observables repeat a name: $names")

    isnothing(priors) || check_priors(priors)

    return CalibrateConfig(
        config_file,
        observables,
        Float64.(pressure_levels),
        n_iterations,
        sample_date_ranges,
        covariance_date_ranges,
        extend,
        spinup,
        output_dir,
        rng_seed,
        priors,
        worker_walltime,
    )

end

"""
    AbstractDataLoader

Supertypes for all data loaders.

`AbstractDataLoader`s have to provide two functions: `available_vars` and
`Base.get`.

`AbstractDataLoader` load preprocessed data as `OutputVar`s that can be broadly
used for any kind of calibration.
"""
abstract type AbstractDataLoader end

"""
    available_vars(data_loader::AbstractDataLoader)

Return the available preprocessed variables in `data_loader`.
"""
available_vars(data_loader::AbstractDataLoader) = data_loader.available_vars

function Base.show(io::IO, data_loader::AbstractDataLoader)
    vars = sort(collect(data_loader.available_vars))
    printstyled(io, nameof(typeof(data_loader)), bold = true, color = :green)
    print(io, ": ")
    print(io, join(vars, ", "))
end

function preprocess(data_loader::AbstractDataLoader, _, ::Val{varname}) where {varname}
    error(
        "No preprocessing function is found for $varname. Add a method for preprocess(::$(typeof(data_loader)), ::Val{:$varname}) that preprocess this variable as a OutputVar",
    )
end

"""
    CompositeDataLoader

A struct for simplifying the process of loading multiple variables from multiple
data loaders.
"""
struct CompositeDataLoader <: AbstractDataLoader
    """A dictionary mapping each variable short name to the `AbstractDataLoader`
    responsible for loading it."""
    varname_to_loaders::Dict{String, AbstractDataLoader}

    """A list of available variables to load."""
    available_vars::Set{String}
end

"""
    CompositeDataLoader(loaders::AbstractDataLoader...; varname_to_loader = Dict())

Construct a `CompositeDataLoader` from multiple data loaders.

The keyword argument `varname_to_loader` is a dictionary mapping variable names
to `AbstractDataLoader`s. When multiple data loaders provide the same variable,
use `varname_to_loader` to specify which loader to use for each variable. If a
variable name is not provided in `varname_to_loader` and multiple loaders
provide the same variable, an error is thrown.

See the example below for how to use `CompositeDataLoader`.

```julia
composite_data_loader = CompositeDataLoader(ERA5DataLoader(), CERESDataLoader())

# If pr is added to ERA5DataLoader, then you can specify to load pr from
# GPCPDataLoader
era5_data_loader = ERA5DataLoader()
gpcp_data_loader = GPCPDataLoader()
composite_data_loader = CompositeDataLoader(
    era5_data_loader,
    gpcp_data_loader;
    varname_to_loader = Dict("pr" => gpcp_data_loader)
)
```
"""
function CompositeDataLoader(loaders::AbstractDataLoader...; varname_to_loader = Dict())
    varname_to_loader = Dict{String, AbstractDataLoader}(varname_to_loader)
    for (varname, loader) in varname_to_loader
        varname in available_vars(loader) ||
            error("$varname is not available in $(nameof(typeof(loader)))")
    end

    # The variable names being the same is a problem when there is ambiguity
    # for which data loaders to get it from
    if length(loaders) != 1
        all_varnames = available_vars.(loaders)
        shared_varnames = intersect(all_varnames...)
        setdiff!(shared_varnames, keys(varname_to_loader))
        if !isempty(shared_varnames)
            error("There are shared variable names between data loaders")
        end
    end

    for loader in loaders
        for varname in available_vars(loader)
            varname in keys(varname_to_loader) && continue
            varname_to_loader[varname] = loader
        end
    end
    return CompositeDataLoader(varname_to_loader, Set(keys(varname_to_loader)))
end

"""
    find_source_loader(loader::CompositeDataLoader, short_name::String)

Find the loader for `short_name` in `loader`.
"""
function find_source_loader(loader::CompositeDataLoader, short_name::String)
    (; varname_to_loaders) = loader
    short_name in keys(varname_to_loaders) || error(
        "$short_name is not available for this composite data loader. Use available_vars to see all available variables",
    )
    return varname_to_loaders[short_name]
end

"""
    get(loader::CompositeDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name`.
"""
function Base.get(loader::CompositeDataLoader, short_name::String)
    (; varname_to_loaders) = loader
    short_name in keys(varname_to_loaders) || error(
        "$short_name is not available to load from this composite data loader. Use available_vars to see all available variables",
    )
    return get(varname_to_loaders[short_name], short_name)
end

function Base.show(io::IO, data_loader::CompositeDataLoader)
    printstyled(io, "CompositeDataLoader", bold = true, color = :green)
    data_loaders = unique(values(data_loader.varname_to_loaders))
    for data_loader in data_loaders
        print(io, "\n  ")
        show(io, data_loader)
    end
end

"""
    ERA5DataLoader

A struct for loading preprocessed ERA5 data as `OutputVar`s.
"""
struct ERA5DataLoader <: AbstractDataLoader
    """A catalog built from NetCDF files, designed to initialize `OutputVar`s.
    See ClimaAnalysis documentation for more information about NCCatalog."""
    catalog::NCCatalog

    """A list of available variables to load."""
    available_vars::Set{String}
end

const ERA5_TO_CLIMA_NAMES =
    ["mslhf" => "hfls", "msshf" => "hfss", "msuwswrf" => "rsus", "msuwlwrf" => "rlus"]
const STANDARD_UNITS = Dict(
    "W m**-2" => "W m^-2",
    "W m-2" => "W m^-2",
    "kg kg**-1" => "unitless",
    "%" => "unitless",
)

"""
    ERA5DataLoader(; era5_to_clima_names = ERA5_TO_CLIMA_NAMES)

Construct a data loader which you can used to load preprocessed ERA5 monthly
time-averaged data in `OutputVar`, where
- the short name, sign of the data, and units match CliMA conventions
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th),
- units match the variables in the output of the CliMA diagnostics.

The ERA5 data comes from the
`era5_monthly_averages_surface_single_level_1979_2024` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/era5_monthly_averages_single_level_1979_2024)
for more information about this artifact.

The keyword argument `era5_to_clima_names` is a vector of pairs mapping
ERA5 name to CliMA name.
"""
function ERA5DataLoader(; era5_to_clima_names = ERA5_TO_CLIMA_NAMES)
    artifact_dir = @clima_artifact"era5_monthly_averages_surface_single_level_1979_2024"
    flux_file = joinpath(
        artifact_dir,
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    )

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, flux_file, era5_to_clima_names...)
    return ERA5DataLoader(catalog, Set(last.(era5_to_clima_names)))
end

"""
    get(loader::ERA5DataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the ERA5
dataset.
"""
function Base.get(loader::ERA5DataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, add it to ERA5_TO_CLIMA_NAMES as a pair mapping ERA5 name to CliMA name and create a new ERA5DataLoader",
    )
    var = get(catalog, short_name; var_kwargs = (shift_by = Dates.firstdayofmonth,))
    return preprocess(loader, var, Val(Symbol(short_name)))
end

"""
    preprocess(::ERA5DataLoader, var, ::Val{varname}) where {varname}

Preprocess `var` with short name `varname`.
"""
preprocess(::ERA5DataLoader, var, ::Val{:hfls}) = _preprocess_var(var; flip_sign = true)
preprocess(::ERA5DataLoader, var, ::Val{:hfss}) = _preprocess_var(var; flip_sign = true)
preprocess(::ERA5DataLoader, var, ::Val{:rsus}) = _preprocess_var(var)
preprocess(::ERA5DataLoader, var, ::Val{:rlus}) = _preprocess_var(var)

function _preprocess_var(var; flip_sign = false)
    short_name = ClimaAnalysis.short_name(var)

    flip_sign && replace!(val -> -val, var)
    issorted(ClimaAnalysis.latitudes(var)) ||
        ClimaAnalysis.reverse_dim!(var, ClimaAnalysis.latitude_name(var))

    # Longitudes in ERA5 dataset is 0 to 360 degrees, but longitudes in CliMA
    # diagnostics range from -180 to 180 degrees
    var = ClimaAnalysis.shift_longitude(var, -180.0, 180.0)

    var_units = ClimaAnalysis.units(var)
    var_units in keys(STANDARD_UNITS) &&
        (var = ClimaAnalysis.set_units(var, STANDARD_UNITS[var_units]))

    # Functions in ClimaAnalysis can mutate the short name. Calibration uses
    # short names to check the observational and simulation data, so we keep the
    # same short name as before.
    ClimaAnalysis.set_short_name!(var, short_name)
    return var
end

"""
    CERESDataLoader

A struct for loading preprocessed CERES data as `OutputVar`s.
"""
struct CERESDataLoader <: AbstractDataLoader
    """A catalog built from NetCDF files, designed to initialize `OutputVar`s.
    See ClimaAnalysis documentation for more information about NCCatalog."""
    catalog::NCCatalog

    """A list of available variables to load."""
    available_vars::Set{String}
end

const CERES_TO_CLIMA_NAMES = [
    "toa_sw_all_mon" => "rsut",
    "toa_sw_clr_t_mon" => "rsutcs",
    "toa_lw_all_mon" => "rlut",
    "toa_lw_clr_t_mon" => "rlutcs",
    "toa_net_all_mon" => "rtmt",
]
const CERES_DERIVED_VARS = Set(["swcre", "lwcre"])

"""
    CERESDataLoader(; ceres_to_clima_names = CERES_TO_CLIMA_NAMES)

Construct a data loader which you can used to load preprocessed CERES monthly
time-averaged TOA radiation data in `OutputVar`, where
- the short name and units match CliMA conventions,
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th).

In addition to the variables in `CERES_TO_CLIMA_NAMES`, the following derived
cloud radiative effect (CRE) variables are available:
- `swcre = rsutcs - rsut` (shortwave cloud radiative effect)
- `lwcre = rlutcs - rlut` (longwave cloud radiative effect)

The CERES data comes from the `radiation_obs` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/radiation_obs)
for more information about this artifact.

The keyword argument `ceres_to_clima_names` is a vector of pairs mapping
CERES name to CliMA name.
"""
function CERESDataLoader(; ceres_to_clima_names = CERES_TO_CLIMA_NAMES)
    artifact_dir = @clima_artifact("radiation_obs")
    ceres_file = joinpath(artifact_dir, "CERES_EBAF_Ed4.2.1_Subset_200003-202512.nc")

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, ceres_file, ceres_to_clima_names...)
    all_vars = Set(last.(ceres_to_clima_names)) ∪ CERES_DERIVED_VARS
    return CERESDataLoader(catalog, all_vars)
end

"""
    get(loader::CERESDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the CERES
dataset.
"""
function Base.get(loader::CERESDataLoader, short_name::String)
    (; available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, add it to CERES_TO_CLIMA_NAMES as a pair mapping CERES name to CliMA name and create a new CERESDataLoader",
    )
    return _get(loader, Val(Symbol(short_name)))
end

# Helper functions to support getting swcre and lwcre
function _get(loader::CERESDataLoader, ::Val{short_name}) where {short_name}
    var = get(
        loader.catalog,
        string(short_name);
        var_kwargs = (shift_by = Dates.firstdayofmonth,),
    )
    return preprocess(loader, var, Val(short_name))
end

function _get(loader::CERESDataLoader, ::Val{:swcre})
    rsutcs = get(loader, "rsutcs")
    rsut = get(loader, "rsut")
    result = rsutcs - rsut
    ClimaAnalysis.set_short_name!(result, "swcre")
    return preprocess(loader, result, Val(:swcre))
end

function _get(loader::CERESDataLoader, ::Val{:lwcre})
    rlutcs = get(loader, "rlutcs")
    rlut = get(loader, "rlut")
    result = rlutcs - rlut
    ClimaAnalysis.set_short_name!(result, "lwcre")
    return preprocess(loader, result, Val(:lwcre))
end

"""
    preprocess(::CERESDataLoader, var, ::Val{varname}) where {varname}

Preprocess `var` with short name `varname`.
"""
preprocess(::CERESDataLoader, var, ::Val{:rsut}) = _preprocess_var(var)
preprocess(::CERESDataLoader, var, ::Val{:rsutcs}) = _preprocess_var(var)
preprocess(::CERESDataLoader, var, ::Val{:rlut}) = _preprocess_var(var)
preprocess(::CERESDataLoader, var, ::Val{:rlutcs}) = _preprocess_var(var)
preprocess(::CERESDataLoader, var, ::Val{:rtmt}) = _preprocess_var(var)
# These variables don't require more postprocessing so we return them as it
preprocess(::CERESDataLoader, var, ::Val{:swcre}) = var
preprocess(::CERESDataLoader, var, ::Val{:lwcre}) = var

"""
    GPCPLoader

A struct for loading preprocessed GPCP data as `OutputVar`s.
"""
struct GPCPDataLoader <: AbstractDataLoader
    catalog::NCCatalog
    available_vars::Set{String}
end

const GPCP_TO_CLIMA_NAMES = ["precip" => "pr"]

"""
    GPCPDataLoader(; gpcp_to_clima_names = GPCP_TO_CLIMA_NAMES)

Construct a data loader which you can used to load preprocessed GPCP monthly
time-averaged precipitation data in `OutputVar`, where
- the short name and sign of the data match CliMA conventions
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th).

The ERA5 data comes from the
`precipitation_obs` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/precipitation_obs)
for more information about this artifact.

Note that the units of precipitation or `pr` are `mm/day` from this dataset. In
CliMA, the units of `pr` are `kg m^-2 s^-1`.

The keyword argument `gpcp_to_clima_names` is a vector of pairs mapping
GPCP name to CliMA name.
"""
function GPCPDataLoader(; gpcp_to_clima_names = GPCP_TO_CLIMA_NAMES)
    artifact_dir = @clima_artifact("precipitation_obs")
    precip_file = joinpath(artifact_dir, "precip.mon.mean.nc")

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, precip_file, gpcp_to_clima_names...)
    return GPCPDataLoader(catalog, Set(last.(gpcp_to_clima_names)))
end

"""
    get(loader::GPCPDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the GPCP
dataset.
"""
function Base.get(loader::GPCPDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, pass it to gpcp_to_clima_names as a pair mapping GPCP name to CliMA name and create a new GPCPDataLoader",
    )
    var = get(catalog, short_name; var_kwargs = (shift_by = Dates.firstdayofmonth,))
    return preprocess(loader, var, Val(Symbol(short_name)))
end

preprocess(::GPCPDataLoader, var, ::Val{:pr}) = _preprocess_var(var, flip_sign = true)

"""
    ERA5PressureLevelDataLoader

A struct for loading preprocessed ERA5 data on pressure levels as `OutputVar`s.
"""
struct ERA5PressureLevelDataLoader <: AbstractDataLoader
    catalog::NCCatalog
    available_vars::Set{String}
end

const ERA5_PRESSURE_LEVEL_TO_CLIMA_NAMES = ["t" => "ta", "q" => "hus", "r" => "hur"]

"""
    ERA5PressureLevelDataLoader(;
        era5_pressure_level_to_clima_names = ERA5_PRESSURE_LEVEL_TO_CLIMA_NAMES)

Construct a data loader which you can used to load preprocessed monthly
time-averaged ERA5 data on pressure levels in `OutputVar`, where
- the short name and sign of the data match CliMA conventions,
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th).

Note that the units of unitless `OutputVar`s are `"unitless"` rather than the
empty string as ClimaAnalysis consider the empty string as missing units.

The ERA5 data comes from the `era5_monthly_averages_pressure_levels_1979_2024`
artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/era5_monthly_averages_pressure_levels_1979_2024)
for more information about this artifact.

The keyword argument `era5_pressure_level_to_clima_names` is a vector of pairs mapping
ERA5 name to CliMA name.
"""
function ERA5PressureLevelDataLoader(;
    era5_pressure_level_to_clima_names = ERA5_PRESSURE_LEVEL_TO_CLIMA_NAMES,
)
    artifact_dir = @clima_artifact"era5_monthly_averages_pressure_levels_1979_2024"
    era5_file =
        joinpath(artifact_dir, "era5_monthly_averages_pressure_levels_197901-202410.nc")

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, era5_file, era5_pressure_level_to_clima_names...)
    return ERA5PressureLevelDataLoader(
        catalog,
        Set(last.(era5_pressure_level_to_clima_names)),
    )
end

"""
    get(loader::ERA5PressureLevelDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the ERA5
pressure levels dataset.
"""
function Base.get(loader::ERA5PressureLevelDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, pass it to era5_pressure_level_to_clima_names as a pair mapping ERA5 name to CliMA name and create a new ERA5PressureLevelDataLoader",
    )
    var = get(catalog, short_name; var_kwargs = (shift_by = Dates.firstdayofmonth,))
    return preprocess(loader, var, Val(Symbol(short_name)))
end

preprocess(::ERA5PressureLevelDataLoader, var, ::Val{:ta}) = _preprocess_var(var)
preprocess(::ERA5PressureLevelDataLoader, var, ::Val{:hus}) = _preprocess_var(var)
function preprocess(::ERA5PressureLevelDataLoader, var, ::Val{:hur})
    var = _preprocess_var(var)
    # Convert from percentages (e.g. 120%) to decimal (1.20)
    var = ClimaAnalysis.convert_units(var, "unitless", conversion_function = x -> 0.01 * x)
    return var
end

"""
    ModisDataLoader

A struct for loading preprocessed MODIS data on pressure levels as `OutputVar`s.
"""
struct ModisDataLoader <: AbstractDataLoader
    catalog::NCCatalog
    available_vars::Set{String}
end

const MODIS_TO_CLIMA_NAMES = ["lwp" => "lwp", "iwp" => "clivi"]

"""
    ModisDataLoader(;
        modis_to_clima_names = MODIS_TO_CLIMA_NAMES,
    )

Construct a data loader which you can used to load preprocessed monthly time-averaged
MODIS data on single levels in `OutputVar`, where
- the short name and sign of the data match CliMA conventions,
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th).

For `lwp` and `clivi`, there are `NaN`s in the data. These `NaN`s was imputted
with mean of the non-`NaN` data for the dataset.

The MODIS data comes from the `modis_lwp_iwp` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/modis_lwp_iwp/)
for more information about this artifact.

The keyword argument `modis_to_clima_names` is a vector of pairs mapping
MODIS name to CliMA name.
"""
function ModisDataLoader(; modis_to_clima_names = MODIS_TO_CLIMA_NAMES)
    artifact_dir = @clima_artifact"modis_lwp_iwp"
    modis_file = joinpath(artifact_dir, "modis_lwp_iwp.nc")

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, modis_file, modis_to_clima_names...)
    return ModisDataLoader(catalog, Set(last.(modis_to_clima_names)))
end

"""
    get(loader::ModisDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the MODIS
dataset.
"""
function Base.get(loader::ModisDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, pass it to modis_to_clima_names as a pair mapping ERA5 name to CliMA name and create a new ModisDataLoader",
    )
    var = get(catalog, short_name; var_kwargs = (shift_by = Dates.firstdayofmonth,))
    return preprocess(loader, var, Val(Symbol(short_name)))
end

function preprocess(::ModisDataLoader, var, ::Union{Val{:clivi}, Val{:lwp}})
    var = _preprocess_var(var)
    not_nans = filter(!isnan, var.data)
    global_mean = sum(not_nans) / length(not_nans)
    @info "$(ClimaAnalysis.short_name(var)): Imputing global mean $global_mean for NaNs"
    replace!(x -> isnan(x) ? global_mean : x, var)
    return var
end

"""
    CalipsoDataLoader

A struct for loading preprocessed CALIPSO/CloudSat cloud fraction data on altitude
levels as `OutputVar`s.
"""
struct CalipsoDataLoader <: AbstractDataLoader
    catalog::NCCatalog
    available_vars::Set{String}
end

const CALIPSO_TO_CLIMA_NAMES = ["cloud_fraction_on_levels" => "cl"]

"""
    CalipsoDataLoader(; calipso_to_clima_names = CALIPSO_TO_CLIMA_NAMES)

Construct a data loader which you can used to load preprocessed monthly CALIPSO/CloudSat
cloud fraction on altitude levels in `OutputVar`, where
- the short name and units match CliMA conventions,
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th),
- values are fractions in [0, 1]

The CALIPSO/CloudSat data comes from the `calipso_cloudsat` artifact, which
repackages the 3S-GEOPROF-COMB dataset. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/calipso_cloudsat)
for more information about this artifact.

The keyword argument `calipso_to_clima_names` is a vector of pairs mapping
CALIPSO/CloudSat variable name to CliMA name.
"""
function CalipsoDataLoader(; calipso_to_clima_names = CALIPSO_TO_CLIMA_NAMES)
    artifact_dir = @clima_artifact("calipso_cloudsat")
    monthly_file = joinpath(artifact_dir, "radarlidar_monthly_2.5x2.5.nc")

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, monthly_file, calipso_to_clima_names...)
    return CalipsoDataLoader(catalog, Set(last.(calipso_to_clima_names)))
end

"""
    get(loader::CalipsoDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the
CALIPSO/CloudSat dataset.
"""
function Base.get(loader::CalipsoDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, pass it to calipso_to_clima_names as a pair mapping CALIPSO/CloudSat name to CliMA name and create a new CalipsoDataLoader",
    )
    var = get(catalog, short_name; var_kwargs = (shift_by = Dates.firstdayofmonth,))
    return preprocess(loader, var, Val(Symbol(short_name)))
end

"""
    preprocess(::CalipsoDataLoader, var, ::Val{varname}) where {varname}

Preprocess `var` with short name `varname`.
"""
function preprocess(::CalipsoDataLoader, var, ::Val{:cl})
    # Select "All cases" doop (first index)
    var = ClimaAnalysis.slice(var; doop = 1, by = ClimaAnalysis.Index())
    var = _preprocess_var(var)
    # Convert from percentage (0–100) to fraction (0–1)
    var = ClimaAnalysis.convert_units(var, "unitless"; conversion_function = x -> x / 100)
    # Impute missing and NaN values with the global mean in one pass
    T = nonmissingtype(eltype(var.data))
    valid = filter(x -> !ismissing(x) && !isnan(x), vec(var.data))
    global_mean = T(sum(valid) / length(valid))
    @info "$(ClimaAnalysis.short_name(var)): Imputing global mean $global_mean for missing/NaN values"
    return ClimaAnalysis.remake(
        var;
        data = map(x -> (ismissing(x) || isnan(x)) ? global_mean : T(x), var.data),
    )
end

"""
    MACDataLoader

A struct for loading preprocessed MAC data on pressure levels as `OutputVar`s.
"""
struct MACDataLoader <: AbstractDataLoader
    catalog::NCCatalog
    available_vars::Set{String}
end

const MAC_TO_CLIMA_NAMES = ["cloudlwp" => "lwp"]

"""
    MACDataLoader(;
        mac_to_clima_names = MAC_TO_CLIMA_NAMES,
    )

Construct a data loader which you can used to load preprocessed monthly time-averaged
MAC data on single levels in `OutputVar`, where
- the short name and sign of the data match CliMA conventions,
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th).

The dataset contains missing values (the file's `_FillValue`). Unlike the
MODIS and CALIPSO loaders, these are left in place and not imputed.

The MAC data comes from the `mac_lwp` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/main/mac_lwp/)
for more information about this artifact.

The keyword argument `mac_to_clima_names` is a vector of pairs mapping
MAC name to CliMA name.
"""
function MACDataLoader(; mac_to_clima_names = MAC_TO_CLIMA_NAMES)
    artifact_dir = @clima_artifact"mac_lwp"
    mac_file = joinpath(artifact_dir, "mac_lwp.nc")

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, mac_file, mac_to_clima_names...)
    return MACDataLoader(catalog, Set(last.(mac_to_clima_names)))
end

"""
    get(loader::MACDataLoader, short_name::String)

Get the preprocessed `OutputVar` with the name `short_name` from the MAC
dataset.
"""
function Base.get(loader::MACDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, pass it to mac_to_clima_names as a pair mapping var name to CliMA name and create a new MACDataLoader",
    )
    var = get(catalog, short_name; var_kwargs = (shift_by = Dates.firstdayofmonth,))
    return preprocess(loader, var, Val(Symbol(short_name)))
end

function preprocess(::MACDataLoader, var, ::Union{Val{:iwp}, Val{:lwp}})
    var = _preprocess_var(var)
    var = ClimaAnalysis.convert_units(var, "kg m^-2", conversion_function = x -> 0.001 * x)
    return drop_empty_times(var)
end

"""
    drop_empty_times(var; max_nan_frac = 0.99)

Drop the times of `var` whose data is more than `max_nan_frac` `NaN`.

A month with essentially no retrievals is not an observation, and it poisons anything that
takes a mask across times: MAC `lwp` has one such month (1988-01, 99.97% `NaN`) against a
median of 49.9% for the rest, which is just its ocean mask. Unifying the `NaN` mask over
the record without dropping it leaves the field entirely empty.

The default is not sensitive. MAC's other months span 47.2% to 56.1% `NaN`, so every
threshold from 0.6 to 0.99 drops that one month and nothing else. Lowering it past the
gap, to 0.55, would drop two further months that are merely below average, which trades
samples for area and is a different decision from dropping empty ones.
"""
function drop_empty_times(var; max_nan_frac = 0.99)
    ClimaAnalysis.has_time(var) || return var
    n_times = length(ClimaAnalysis.times(var))
    keep = filter(1:n_times) do index
        slice = ClimaAnalysis.view_select(var; by = ClimaAnalysis.Index(), time = index)
        count(isnan, slice.data) / length(slice.data) <= max_nan_frac
    end
    length(keep) == n_times && return var
    isempty(keep) && error(
        "Every time of $(ClimaAnalysis.short_name(var)) is more than $(100max_nan_frac)% NaN",
    )
    @info "Dropping $(n_times - length(keep)) empty time(s) from $(ClimaAnalysis.short_name(var))"
    return ClimaAnalysis.select(var; by = ClimaAnalysis.Index(), time = keep)
end

"""
    default_data_loader()

One loader for every observation a calibration can use, by short name: ERA5 pressure-level
`ta`, `hus`, `hur`; ERA5 surface `hfls`, `hfss`, `rsus`, `rlus`; CERES `rsut`, `rlut`, their
clear-sky versions, `swcre`, `lwcre` and `rtmt`; MODIS `clivi`; MAC `lwp`; GPCP `pr`.
To add an observation, add its loader here. Each field needs a model counterpart in
`model_var`.
"""
function default_data_loader()
    mac = MACDataLoader()
    return CompositeDataLoader(
        ERA5PressureLevelDataLoader(),
        ERA5DataLoader(),
        CERESDataLoader(),
        ModisDataLoader(),
        mac,
        GPCPDataLoader();
        # MODIS also provides `lwp`; MAC's is the one used.
        varname_to_loader = Dict("lwp" => mac),
    )
end

"""
    model_var(simdir, short_name)

The model's monthly means of `short_name` from the diagnostics in `simdir`, as a vector of
`OutputVar`s, in the units and sign convention of the observation of the same name in
`default_data_loader`. A field output on both pressure and altitude coordinates gives one
`OutputVar` for each.

Fields the model does not output directly, or outputs in other units, have their own method
on `Val(Symbol(short_name))`: add one alongside the observation's loader for a new field.
"""
model_var(simdir, short_name::String) = model_var(simdir, Val(Symbol(short_name)))

function model_var(simdir, ::Val{short_name}) where {short_name}
    name = string(short_name)
    coord_types = ClimaAnalysis.available_coord_types(
        simdir;
        short_name = name,
        reduction = "average",
        period = "1M",
    )
    return [
        get(simdir; short_name = name, reduction = "average", period = "1M", coord_type) for
        coord_type in coord_types
    ]
end

monthly_mean(simdir, name) =
    get(simdir; short_name = name, reduction = "average", period = "1M")

function toa_flux(var, short_name)
    ClimaAnalysis.set_short_name!(var, short_name)
    return [ClimaAnalysis.set_units(var, "W m^-2")]
end

model_var(simdir, ::Val{:swcre}) =
    toa_flux(monthly_mean(simdir, "rsutcs") - monthly_mean(simdir, "rsut"), "swcre")
model_var(simdir, ::Val{:lwcre}) =
    toa_flux(monthly_mean(simdir, "rlutcs") - monthly_mean(simdir, "rlut"), "lwcre")
# Net downward radiation at the top of the atmosphere.
model_var(simdir, ::Val{:rtmt}) = toa_flux(
    monthly_mean(simdir, "rsdt") - monthly_mean(simdir, "rsut") -
    monthly_mean(simdir, "rlut"),
    "rtmt",
)

# The model's precipitation is kg m^-2 s^-1. GPCP is mm/day, and its loader gives it the
# model's sign convention, so only the magnitude is converted.
function model_var(simdir, ::Val{:pr})
    pr = monthly_mean(simdir, "pr")
    pr = ClimaAnalysis.remake(pr; data = pr.data .* 86400)
    return [ClimaAnalysis.set_units(pr, "mm/day")]
end

"""
    regrid_to_model_levels(obs_var, model_altitudes)

Regrid `obs_var` from dataset altitude levels to `model_altitudes` using a
Gaussian-weighted average along the altitude dimension.

For each model level `k` with altitude `z_k`, the regridded value is:

    obs(z_k) = sum_i w_i(k) * obs(z_i)

where the sum is over all dataset levels `i` with altitude `z_i`, and the
normalized weight is:

    w_i(k) = exp(-(z_i - z_k)^2 / (2 * sigma_k^2)) / Z_k
    Z_k     = sum_j exp(-(z_j - z_k)^2 / (2 * sigma_k^2))

The bandwidth `sigma_k` is half the spacing between neighboring model levels:

    sigma_k = 0.5 * (z_{k+1} - z_{k-1})   (interior levels)
    sigma_1 = 0.5 * (z_2 - z_1)            (lower boundary)
    sigma_K = 0.5 * (z_K - z_{K-1})        (upper boundary)

All other dimensions (lon, lat, time, …) are left unchanged; only the altitude
axis is replaced by `model_altitudes`.

`obs_var` must have an altitude dimension. `model_altitudes` should be in the
same units as the altitude coordinate of `obs_var` (meters above mean sea
level for the CALIPSO/CloudSat dataset).
"""
function regrid_to_model_levels(obs_var::ClimaAnalysis.OutputVar, model_altitudes)
    ClimaAnalysis.has_altitude(obs_var) || error(
        "obs_var with short name $(short_name(obs_var)) is missing the altitude dimension",
    )
    alt_name = ClimaAnalysis.altitude_name(obs_var)
    alt_axis = obs_var.dim2index[alt_name]
    z_i = ClimaAnalysis.altitudes(obs_var)
    z_k = model_altitudes
    K = length(z_k)

    # Compute sigma for each model level (half the central finite difference)
    sigmas = map(1:K) do k
        if k == 1
            0.5 * (z_k[2] - z_k[1])
        elseif k == K
            0.5 * (z_k[K] - z_k[K - 1])
        else
            0.5 * (z_k[k + 1] - z_k[k - 1])
        end
    end

    # Allocate output array with model altitude size along alt_axis
    out_shape = ntuple(i -> i == alt_axis ? K : size(obs_var.data, i), ndims(obs_var.data))
    out_data = zeros(eltype(obs_var.data), out_shape)

    # For each model level k, compute normalized Gaussian weights over all
    # dataset levels and accumulate the weighted sum into the corresponding
    # slice of out_data
    for k in 1:K
        weights = @. exp(-(z_i - z_k[k])^2 / (2 * sigmas[k]^2))
        weights ./= sum(weights)
        out_k = selectdim(out_data, alt_axis, k)
        for (i, w) in enumerate(weights)
            out_k .+= w .* selectdim(obs_var.data, alt_axis, i)
        end
    end

    new_dims = copy(obs_var.dims)
    new_dims[alt_name] = model_altitudes
    return ClimaAnalysis.remake(obs_var; dims = new_dims, data = out_data)
end

"""
    update_timespan!(
        config_dict,
        start_date::Dates.DateTime,
        end_date::Dates.DateTime
    )

Update "start\\_date" and "t\\_end" in `config_dict` to match `start_date` and
`end_date`.

The `start\\_date` and `end\\_date` are converted to strings and the keys
"start_date" and "t_end" in `config_dict` are updated accordingly.
"""
function update_timespan!(config_dict, start_date::Dates.DateTime, end_date::Dates.DateTime)
    start_date_str = Utilities.format_start_date(start_date)
    config_dict["start_date"] = start_date_str
    sim_length = Dates.Second(end_date - start_date)
    config_dict["t_end"] = "$(sim_length.value)secs"
    return nothing
end

"""
    add_parameter_filepath!(config_dict, parameter_filepath)

Add the parameter toml file at `parameter_filepath` to `config_dict`.
"""
function add_parameter_filepath!(config_dict, parameter_filepath)
    parameter_filepaths = get!(config_dict, "coupler_toml", String[])
    push!(parameter_filepaths, parameter_filepath)
    return nothing
end

"""
    setup_and_emulate_diagnostics(config_dict::AbstractDict; dt = "86400secs")

Set up and emulate producing the diagnostics of the coupled model simulation
specified by the `config_dict`. Returns the `CoupledSimulation` after simulating
the diagnostics.

This is analogous to [`SimCoordinator.setup_and_run`](@ref) and is meant a
drop-in replacement for starting a simulation that will produce only
diagnostics. Note this function is primarily intended for testing, not running
or viewing meaningful simulations.

By default, the size of the time step will be set to one day.

!!! note "Atmos diagnostics"
    This only supports simulating the creation of atmos diagnostics.

!!! note "Emulation"
    Since no actual time steps are taken, the values in the `CoupledSimulation`
    might be nonsensical.
"""
function setup_and_emulate_diagnostics(config_dict::AbstractDict; dt = "86400secs")
    config_dict["dt"] = dt
    cs = ClimaCoupler.Interfacer.CoupledSimulation(config_dict)
    emulate_diagnostics!(cs)
    return cs
end

"""
    emulate_diagnostics!(cs::ClimaCoupler.Interfacer.CoupledSimulation)

Emulate producing diagnostics for `cs`.

This only produces diagnostics for the atmos model because the calibration of
the coupled model currently uses only atmos variables.
"""
function emulate_diagnostics!(cs::ClimaCoupler.Interfacer.CoupledSimulation)
    t_end = cs.tspan[end]
    emulate_diagnostics!(cs.model_sims.atmos_sim, t_end)
end

"""
    emulate_diagnostics!(sim, t_end)

Emulate producing diagnostics for the component model `sim`.

There is only support for emulating diagnostics of the `ClimaAtmosSimulation`.
"""
function emulate_diagnostics!(sim, t_end)
    (; integrator, output_writers) = sim
    diagnostics_handler =
        integrator.callback.discrete_callbacks[end].affect!.diagnostics_handler
    while sim.integrator.t < t_end
        sim.integrator.t += integrator.dt
        ClimaDiagnostics.orchestrate_diagnostics(integrator, diagnostics_handler)
    end
    # Since the simulation never actually ran, we need to manually close the
    # writers ourselves
    foreach(close, output_writers)
    return nothing
end

end
