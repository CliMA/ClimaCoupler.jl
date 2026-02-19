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
