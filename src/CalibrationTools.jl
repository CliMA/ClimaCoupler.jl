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

function _show(io::IO, data_loader::AbstractDataLoader, data_loader_name::String)
    vars = sort(collect(data_loader.available_vars))
    printstyled(io, data_loader_name, bold = true, color = :green)
    print(io, ": ")
    print(io, join(vars, ", "))
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
const STANDARD_UNITS =
    Dict("W m**-2" => "W m^-2", "kg kg**-1" => "unitless", "%" => "unitless")

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
    get(loader::ERA5DataLoader, short_name)

Get the preprocessed `OutputVar` with the name `short_name` from the ERA5
dataset.
"""
function Base.get(loader::ERA5DataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, add it to ERA5_TO_CLIMA_NAMES as a pair mapping ERA5 name to CliMA name and create a new ERA5DataLoader",
    )
    var = ClimaAnalysis.Catalog.get(
        catalog,
        short_name;
        var_kwargs = (shift_by = Dates.firstdayofmonth,),
    )
    return preprocess(loader, var, Val(Symbol(short_name)))
end

"""
    preprocess(::ERA5DataLoader, var, ::Val{varname}) where {varname}

Preprocess `var` with short name `varname`.
"""
function preprocess(::ERA5DataLoader, _, ::Val{varname}) where {varname}
    error(
        "No preprocessing function is found for $varname. Add a method for preprocess(var, ::Val{:$varname}) that preprocess this variable as a OutputVar",
    )
end

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

Base.show(io::IO, data_loader::ERA5DataLoader) = _show(io, data_loader, "ERA5DataLoader")

"""
    GPCPLoader

A struct for loading preprocessed GPCP data as `OutputVar`s.
"""
struct GPCPDataLoader <: AbstractDataLoader
    catalog::NCCatalog
    available_vars::Set{String}
end

"""
    GPCPDataLoader(; gpcp_to_clima_names = ["precip" => "pr"])

Construct a data loader which you can load preprocessed GPCP monthly
time-averaged precipitation data in `OutputVar`, where
- the short name and sign of the data match CliMA conventions
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th).

The ERA5 data comes from the
`precipitation_obs` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/f09429b0a149fd9f780bfc1833e4675c9fab573e/precipitation_obs)
for more information about this artifact.

The keyword argument `era5_to_clima_names` is a vector of pairs mapping
ERA5 name to CliMA name.
"""
function GPCPDataLoader(; gpcp_to_clima_names = ["precip" => "pr"])
    artifact_dir = @clima_artifact("precipitation_obs")
    precip_file = joinpath(artifact_dir, "precip.mon.mean.nc")

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, precip_file, gpcp_to_clima_names...)
    return GPCPDataLoader(catalog, Set(last.(gpcp_to_clima_names)))
end

"""
    get(loader::GPCPDataLoader, short_name)

Get the preprocessed `OutputVar` with the name `short_name` from the GPCP
dataset.

Note that the units of precipitation or `pr` are `mm/day` from this dataset. In
CliMA, the units of `pr` are `kg m^-2 s^-1`.
"""
function Base.get(loader::GPCPDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, pass it to gpcp_to_clima_names as a pair mapping GPCP name to CliMA name and create a new GPCPDataLoader",
    )
    var = ClimaAnalysis.Catalog.get(
        catalog,
        short_name;
        var_kwargs = (shift_by = Dates.firstdayofmonth,),
    )
    return preprocess(loader, var, Val(Symbol(short_name)))
end

preprocess(::GPCPDataLoader, var, ::Val{:pr}) = _preprocess_var(var, flip_sign = true)

Base.show(io::IO, data_loader::GPCPDataLoader) = _show(io, data_loader, "GPCPDataLoader")

"""
    ERA5PressureLevelDataLoader

A struct for loading preprocessed ERA5 data on pressure levels as `OutputVar`s.
"""
struct ERA5PressureLevelDataLoader <: AbstractDataLoader
    catalog::NCCatalog
    available_vars::Set{String}
end

ERA5_PRESSURE_LEVEL_TO_CLIMA_NAMES = ["t" => "ta", "q" => "hus", "r" => "hur"]

"""
    ERA5PressureLevelDataLoader(;
        era5_pressure_level_to_clima_names = ERA5_PRESSURE_LEVEL_TO_CLIMA_NAMES)

Construct a data loader which you can load preprocessed monthly
time-averaged ERA5 data on pressure levels in `OutputVar`, where
- the short name and sign of the data match CliMA conventions,
- the latitudes are shifted to be -180 to 180 degrees,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th).

The ERA5 data comes from the
`era5_monthly_averages_pressure_levels_1979_2024` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts/tree/f09429b0a149fd9f780bfc1833e4675c9fab573e/era5_monthly_averages_pressure_levels_1979_2024)
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

Note that the units of unitless `OutputVar`s are `"unitless"` rather than the
empty string as ClimaAnalysis consider the empty string as missing units.
"""
function Base.get(loader::ERA5PressureLevelDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, pass it to era5_pressure_level_to_clima_names as a pair mapping ERA5 name to CliMA name and create a new ERA5PressureLevelDataLoader",
    )
    var = ClimaAnalysis.Catalog.get(
        catalog,
        short_name;
        var_kwargs = (shift_by = Dates.firstdayofmonth,),
    )
    return preprocess(loader, var, Val(Symbol(short_name)))
end

preprocess(::ERA5PressureLevelDataLoader, var, ::Val{:ta}) = _preprocess_var(var)
preprocess(::ERA5PressureLevelDataLoader, var, ::Val{:hus}) = _preprocess_var(var)
function preprocess(::ERA5PressureLevelDataLoader, var, ::Val{:hur})
    var = _preprocess_var(var)
    # Convert from percentages (e.g. 90%) to decimal (0.90)
    var = ClimaAnalysis.convert_units(var, "unitless", conversion_function = x -> 0.01 * x)
    return var
end

Base.show(io::IO, data_loader::ERA5PressureLevelDataLoader) =
    _show(io, data_loader, "ERA5PressureLevelDataLoader")

"""
    CERESDataLoader

A struct for loading preprocessed CERES EBAF radiation data as `OutputVar`s.
"""
struct CERESDataLoader <: AbstractDataLoader
    """A catalog built from NetCDF files, designed to initialize `OutputVar`s.
    See ClimaAnalysis documentation for more information about NCCatalog."""
    catalog::NCCatalog

    """A list of available variables to load."""
    available_vars::Set{String}
end

const CERES_TO_CLIMA_NAMES = [
    "solar_mon" => "rsdt",
    "toa_sw_all_mon" => "rsut",
    "toa_lw_all_mon" => "rlut",
    "toa_sw_clr_t_mon" => "rsutcs",
    "toa_lw_clr_t_mon" => "rlutcs",
    "sfc_sw_down_all_mon" => "rsds",
    "sfc_sw_up_all_mon" => "rsus",
    "sfc_lw_down_all_mon" => "rlds",
    "sfc_lw_up_all_mon" => "rlus",
    "sfc_sw_down_clr_t_mon" => "rsdscs",
    "sfc_sw_up_clr_t_mon" => "rsuscs",
    "sfc_lw_down_clr_t_mon" => "rldscs",
]
const CERES_TO_CLIMA_UNITS = Dict("W m-2" => "W m^-2")

"""
    CERESDataLoader(; ceres_to_clima_names = CERES_TO_CLIMA_NAMES)

Construct a data loader which you can load preprocessed CERES EBAF monthly
time-averaged radiation data in `OutputVar`, where
- the short name and units match CliMA conventions,
- the times are at the start of the time period (e.g. the time average of
  January is on the first of January instead of January 15th),
- units match the variables in the output of the CliMA diagnostics.

The CERES data comes from the `radiation_obs` artifact. See
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts)
for more information about this artifact.

The keyword argument `ceres_to_clima_names` is a vector of pairs mapping
CERES name to CliMA name.
"""
function CERESDataLoader(; ceres_to_clima_names = CERES_TO_CLIMA_NAMES)
    artifact_dir = @clima_artifact"radiation_obs"
    radiation_file = joinpath(
        artifact_dir,
        "CERES_EBAF_Ed4.2_Subset_200003-201910.nc",
    )

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(catalog, radiation_file, ceres_to_clima_names...)
    return CERESDataLoader(catalog, Set(last.(ceres_to_clima_names)))
end

"""
    available_vars(data_loader::CERESDataLoader)

Return the available preprocessed variables in `data_loader`.
"""
available_vars(data_loader::CERESDataLoader) = data_loader.available_vars

"""
    get(loader::CERESDataLoader, short_name)

Get the preprocessed `OutputVar` with the name `short_name` from the CERES
dataset.
"""
function Base.get(loader::CERESDataLoader, short_name::String)
    (; catalog, available_vars) = loader
    short_name in available_vars || error(
        "$short_name is not available to load. To add this variable, add it to CERES_TO_CLIMA_NAMES as a pair mapping CERES name to CliMA name and create a new CERESDataLoader",
    )
    var = ClimaAnalysis.Catalog.get(
        catalog,
        short_name;
        var_kwargs = (shift_by = Dates.firstdayofmonth,),
    )
    return preprocess(loader, var, Val(Symbol(short_name)))
end

for short_name in (:rsdt, :rsut, :rlut, :rsutcs, :rlutcs, :rsds, :rsus, :rlds, :rlus, :rsdscs, :rsuscs, :rldscs)
    @eval preprocess(::CERESDataLoader, var, ::Val{$(QuoteNode(short_name))}) =
        _preprocess_var(var)
end

Base.show(io::IO, data_loader::CERESDataLoader) =
    _show(io, data_loader, "CERESDataLoader")

end
