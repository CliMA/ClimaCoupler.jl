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
    ERA5DataLoader

A struct for loading preprocessed ERA5 data as `OutputVar`s.
"""
struct ERA5DataLoader
    """A catalog built from NetCDF files, designed to initialize `OutputVar`s.
    See ClimaAnalysis documentation for more information about NCCatalog."""
    catalog::NCCatalog

    """A list of available variables to load."""
    available_vars::Set{String}
end

const ERA5_TO_CLIMA_NAMES =
    ["mslhf" => "hfls", "msshf" => "hfss", "msuwswrf" => "rsus", "msuwlwrf" => "rlus"]
const ERA5_TO_CLIMA_UNITS = Dict("W m**-2" => "W m^-2")

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
    available_vars(data_loader::ERA5DataLoader)

Return the available preprocessed variables in `data_loader`.
"""
available_vars(data_loader::ERA5DataLoader) = data_loader.available_vars

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
    var_units in keys(ERA5_TO_CLIMA_UNITS) &&
        (var = ClimaAnalysis.set_units(var, ERA5_TO_CLIMA_UNITS[var_units]))

    # Functions in ClimaAnalysis can mutate the short name. Calibration uses
    # short names to check the observational and simulation data, so we keep the
    # same short name as before.
    ClimaAnalysis.set_short_name!(var, short_name)
    return var
end

function Base.show(io::IO, data_loader::ERA5DataLoader)
    vars = sort(collect(data_loader.available_vars))
    printstyled(io, "ERA5DataLoader", bold = true, color = :green)
    print(io, ": ")
    print(io, join(vars, ", "))
end

end
