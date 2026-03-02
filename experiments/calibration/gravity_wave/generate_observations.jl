import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCoupler
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2
import ClimaDiagnostics
import ClimaCore
import ClimaUtilities.ClimaArtifacts.@clima_artifact
import Pkg.Artifacts
# Access CalibrateConfig
include(joinpath(@__DIR__, "run_calibration.jl"))
include(joinpath(@__DIR__, "observation_utils.jl"))

# Path to the ClimaEarth Artifacts.toml file
const CLIMAEARTH_ARTIFACTS_TOML =
    joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "Artifacts.toml")

"""
    load_var(filepath, short_name; varname=nothing, flip_sign=false, transform_dates=false)

Helper function to load and process an ERA5 variable with common transformations.
"""
function load_var(
    filepath,
    short_name;
    varname = nothing,
    flip_sign = false,
    transform_dates = nothing,
)
    # TODO: Don't assume monthly data
    var =
        isnothing(varname) ? OutputVar(filepath) :
        OutputVar(filepath, varname; shift_by = Dates.firstdayofmonth)

    flip_sign && (var = -var)

    var.attributes["short_name"] = short_name

    if !isnothing(transform_dates)
        var = ClimaAnalysis.Var._shift_by(var, date -> date - transform_dates)
    end
    if !issorted(latitudes(var))
        var = reverse_dim(var, latitude_name(var))
    end

    return var
end

"""
    load_vars(obsdir, short_names)

Load NetCDF files belonging to `short_names` in `obsdir` as `OutputVar`s.
"""
function load_vars()
    # Use Artifacts API directly with explicit path to ClimaEarth Artifacts.toml
    # This is necessary because @clima_artifact searches upward from the source
    # file, which won't find experiments/ClimaEarth/Artifacts.toml

    artifact_path = Artifacts.ensure_artifact_installed(
        "era5_monthly_averages_surface_single_level_1979_2024",
        CLIMAEARTH_ARTIFACTS_TOML,
    )
    flux_file = joinpath(
        artifact_path,
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    )

    lhf = load_var(flux_file, "hfls"; varname = "mslhf", flip_sign = true)
    shf = load_var(flux_file, "hfss"; varname = "msshf", flip_sign = true)
    rsus = load_var(flux_file, "rsus"; varname = "msuwswrf")
    rlus = load_var(flux_file, "rlus"; varname = "msuwlwrf")

    return [lhf, shf, rsus, rlus]
end

"""
    load_pfull_vars()

Load ERA5 pressure-level variables (temperature, zonal wind, meridional wind)
for 3D calibration observations.
"""
function load_pfull_vars()
    artifact_path = Artifacts.ensure_artifact_installed(
        "era5_monthly_averages_pressure_levels_1979_2024",
        CLIMAEARTH_ARTIFACTS_TOML,
    )
    pfull_file = joinpath(
        artifact_path,
        "era5_monthly_averages_pressure_levels_197901-202410.nc",
    )

    # Load with model short_name, ERA5 varname
    ta = load_var(pfull_file, "ta"; varname = "t")
    ua = load_var(pfull_file, "ua"; varname = "u")
    va = load_var(pfull_file, "va"; varname = "v")

    # Convert pressure dimension from Pa to hPa
    vars = map([ta, ua, va]) do var
        ClimaAnalysis.Var.convert_dim_units(
            var,
            "pressure_level",
            "hPa";
            conversion_function = x -> 0.01 * x,
        )
    end

    return vars
end

"""
    preprocess_vars(vars)

Preprocess each OutputVar in `vars` by resampling to the model grid
and setting units.

Automatically detects 3D pressure-level variables and applies appropriate resampling.
"""
function preprocess_vars(vars)
    resample_2d = resampled_lonlat(CALIBRATE_CONFIG.config_file; include_pressure = false)
    resample_3d = resampled_lonlat(CALIBRATE_CONFIG.config_file; include_pressure = true)

    vars = map(vars) do var
        # Check if variable has a pressure dimension
        if ClimaAnalysis.has_pressure(var)
            var = resample_3d(var)
        else
            var = resample_2d(var)
        end
        var = set_units(var, var_units[short_name(var)])
    end

    return vars
end

function make_observation_vector(vars, sample_date_ranges)
    covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
        scalar = 5.0,
        use_latitude_weights = true,
    )
    obs_vec = map(sample_date_ranges) do sample_date_range
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            vars,
            first(sample_date_range),
            last(sample_date_range),
        )
    end
    return obs_vec
end

"""
    resampled_lonlat(config_file; include_pressure = false)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.

If `include_pressure` is true, also resample on standard ERA5 pressure levels.
"""
function resampled_lonlat(config_file; include_pressure = false)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(CALIBRATE_CONFIG.config_file)
    if !isnothing(config_dict["netcdf_interpolation_num_points"])
        (nlon, nlat, nlev) = tuple(config_dict["netcdf_interpolation_num_points"]...)
    else
        cs = CoupledSimulation(config_file)
        center_space = cs.model_sims.atmos_sim.domain.center_space
        (nlon, nlat, nlev) = ClimaDiagnostics.Writers.default_num_points(center_space)
        stretch = center_space.grid.vertical_grid.topology.mesh.stretch
        dz_bottom = center_space.grid.vertical_grid.topology.mesh.faces[2].z
        z = range(dz_bottom, ClimaCore.Spaces.z_max(center_space), nlev)
    end
    lon = range(-180, 180, nlon)
    lat = range(-90, 90, nlat)

    if include_pressure
        # Standard ERA5 pressure levels in hPa
        # Note: 1000 and 925 hPa excluded because model z-windowing removes near-surface levels
        pressure_level = [850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 250.0, 200.0]
        return var -> resampled_as(var; lon, lat, pressure_level)
    else
        return var -> resampled_as(var; lon, lat)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file
    @info "Generating observations for $short_names"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    # Load 2D surface variables
    unprocessed_vars = load_vars()

    # Load 3D pressure-level variables
    unprocessed_pfull_vars = load_pfull_vars()

    # Combine all variables
    all_unprocessed_vars = vcat(unprocessed_vars, unprocessed_pfull_vars)

    # Filter to only requested short_names
    all_unprocessed_vars = filter(
        var -> ClimaAnalysis.short_name(var) in short_names,
        all_unprocessed_vars,
    )

    preprocessed_vars = preprocess_vars(all_unprocessed_vars)

    JLD2.save_object(
        joinpath(
            pkgdir(ClimaCoupler),
            "experiments/calibration/gravity_wave/preprocessed_vars.jld2",
        ),
        preprocessed_vars,
    )
    observation_vector = make_observation_vector(preprocessed_vars, sample_date_ranges)
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/gravity_wave/obs_vec.jld2"),
        observation_vector,
    )
end
