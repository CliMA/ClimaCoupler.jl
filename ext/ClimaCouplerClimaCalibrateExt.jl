module ClimaCouplerClimaCalibrateExt

import ClimaCoupler: Calibrate
import ClimaAnalysis: NCCatalog

# TODO: Define a struct for loading in data

struct ERA5DataLoader
    catalog::NCCatalog
end

"""
    ERA5DataLoader

Construct a ERA5 data loader which you can load `OutputVar` from.
"""
function Calibrate.ERA5DataLoader()
    artifact_path = Artifacts.ensure_artifact_installed(
        "era5_monthly_averages_surface_single_level_1979_2024",
        CLIMAEARTH_ARTIFACTS_TOML,
    )
    flux_file = joinpath(
        artifact_path,
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    )

    catalog = NCCatalog()
    ClimaAnalysis.add_file!(
        catalog,
        flux_file,
        "mslhf" => "hfls",
        "msshf" => "hfss",
        "msuwswrf" => "rsus",
        "msuwlwrf" => "rlus",
    )
    ERA5DataLoader(catalog)
end

"""
    load_var(loader::ERA5DataLoader, short_name)

Load
"""
function Calibrate.load_var(loader::ERA5DataLoader, short_name)
    (; catalog) = loader
    var = get(catalog, short_name, var_kwargs = (shift_by = Dates.firstdayofmonth,))

end

# TODO: Finish these
function _preprocess(var, ::hfls)
    var = -var
end

function _preprocess(var, ::hfss)
    var = -var
end

function _preprocess(var, ::rsus)
end

function _preprocess(var, ::rlus)
end

function _preprocess(var)
    if !issorted(latitudes(var))
        var = reverse_dim(var, latitude_name(var))
    end
    return var
end

end
