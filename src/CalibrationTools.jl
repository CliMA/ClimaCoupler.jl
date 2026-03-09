"""
    CalibrationTools

This module contains utilities to simplify calibration experiments of a coupled
simulation.
"""
module CalibrationTools

import Dates

import ClimaAnalysis
import ClimaAnalysis: NCCatalog
import ClimaUtilities.ClimaArtifacts: @clima_artifact

"""
    struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period}

A configuration struct for keeping track of multiple fields that are of interest
to a user running calibration, or that are needed in multiple places (e.g., for
ensemble members and generating observations).
"""
struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period}
    "Configuration file to use for ClimaCoupler simulation"
    config_file::String

    "The short names of the observations used for calibration. The short names
    should match the same names used for the diagnostics."
    short_names::Vector{String}

    "The size of the minibatch for each iteration"
    minibatch_size::Int64

    "The number of iterations to run the calibration for"
    n_iterations::Int64

    "The date ranges of the samples for calibration and used to determine the
    start and end dates of a simulation for each iteration of calibration"
    sample_date_ranges::Vector{NTuple{2, Dates.DateTime}}

    "The amount of time to run a simulation after the last date of the
    minibatch"
    extend::EXTEND

    "The amount of time to run a simulation before the first date of the
    minibatch"
    spinup::SPINUP

    "The directory to store the iterations and members of the calibration."
    output_dir::String

    "An integer value for ensuring calibrations are the same between multiple
    calibrations with the same settings"
    rng_seed::Int64
end

"""
    CalibrateConfig(;
        config_file,
        short_names::Vector{String},
        minibatch_size::Integer,
        n_iterations::Integer,
        sample_date_ranges,
        extend::Dates.Period,
        spinup::Dates.Period,
        output_dir,
        rng_seed = 42,
    )

Initializes a `CalibrateConfig` which contains values needed in multiple places
during calibration.

Keyword arguments
=====================

- `config_file`: Configuration file to use for ClimaCoupler simulation.

- `short_names`: Short names of the observations.

- `minibatch_size`: The size of the minibatch for each iteration.

- `n_iterations`: The number of iterations to run the calibration for.

- `sample_date_ranges`: The date ranges for each sample. The dates should be the
  same as found in the time series data of the observations.

- `extend`: The amount of time to run the simulation after the end date
  determined by `sample_date_ranges`. For seasonal averages, `extend` should be
  `Dates.Month(3)` and for monthly averages, `extend` should be
  `Dates.Month(1)`.

- `spinup`: The amount of time to run the simulation before the start date
  determined by `sample_date_ranges`.

- `output_dir`: The location to save the calibration at.

- `rng_seed`: An integer to ensure that calibration runs with the same settings
  are the same.
"""
function CalibrateConfig(;
    config_file,
    short_names::Vector{String},
    minibatch_size::Integer,
    n_iterations::Integer,
    sample_date_ranges,
    extend::Dates.Period,
    spinup::Dates.Period,
    output_dir,
    rng_seed = 42,
)
    isfile(config_file) || error("Configuration file ($config_file) does not exist")
    isempty(short_names) && error("Cannot run calibration with no observation short names")
    isempty(sample_date_ranges) &&
        error("Cannot run calibration with no date ranges for the samples")

    sign(extend) == -1 && error("The period to extend ($extend) by should be positive")
    sign(spinup) == -1 && error("The period to spin up ($spinup) by should be positive")

    sample_date_ranges = [
        (Dates.DateTime(date_pair[1]), Dates.DateTime(date_pair[2])) for
        date_pair in sample_date_ranges
    ]

    for (start_date, stop_date) in sample_date_ranges
        start_date <= stop_date || error(
            "The start date ($start_date) should be before the stop date ($stop_date)",
        )
    end

    minibatch_size > 0 || error("The minibatch size ($minibatch_size) should be positive")
    n_iterations > 0 || error("The number of iterations ($n_iterations) should be positive")

    num_samples = length(sample_date_ranges)
    minibatch_size > num_samples && error(
        "The minibatch size is $minibatch_size, but the number of samples is $num_samples",
    )

    remaining = num_samples % minibatch_size
    remaining == 0 || @warn(
        "Number of samples is not divisible by the minibatch size; the last $remaining samples may be missing when running the calibration"
    )

    return CalibrateConfig(
        config_file,
        short_names,
        minibatch_size,
        n_iterations,
        sample_date_ranges,
        extend,
        spinup,
        output_dir,
        rng_seed,
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

Construct a data loader which you can load preprocessed ERA5 monthly
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
]
const CERES_DERIVED_VARS = Set(["swcre", "lwcre"])

"""
    CERESDataLoader(; ceres_to_clima_names = CERES_TO_CLIMA_NAMES)

Construct a data loader which you can load preprocessed CERES monthly
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
    ceres_file = joinpath(artifact_dir, "CERES_EBAF_Ed4.2_Subset_200003-201910.nc")

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

Construct a data loader which you can load preprocessed GPCP monthly
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

Construct a data loader which you can load preprocessed monthly
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

Construct a data loader which you can load preprocessed monthly time-averaged
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
    @info "$(ClimaAnalysis.short_name(var)): Imputting $global_mean for NaNs"
    replace!(x -> isnan(x) ? global_mean : x, var)
    return var
end

"""
    update_timespan!(
        config_dict,
        start_date::Dates.DateTime,
        end_date::Dates.DateTime
    )

Update "start_date" and "t_end" in `config_dict` to match `start_date` and
`end_date`.

The `start_date` and `end_date` are converted to strings and the keys
"start_date" and "t_end" in `config_dict` are updated accordingly. Note that any
precision beyond days (e.g. hours, seconds, etc.) are not used for setting the
start date.
"""
function update_timespan!(config_dict, start_date::Dates.DateTime, end_date::Dates.DateTime)
    start_date_str = Dates.format(start_date, "yyyymmdd")
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

end
