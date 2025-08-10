import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2
import ClimaDiagnostics
import ClimaCore
import ClimaUtilities.ClimaArtifacts.@clima_artifact
# Access CalibrateConfig
include(joinpath(@__DIR__, "run_calibration.jl"))
include(joinpath(@__DIR__, "observation_utils.jl"))

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
    flux_file = joinpath(
        @clima_artifact("era5_monthly_averages_surface_single_level_1979_2024"),
        "era5_monthly_averages_surface_single_level_197901-202410.nc",
    )

    lhf = load_var(flux_file, "hfls"; varname = "mslhf")
    shf = load_var(flux_file, "hfss"; varname = "msshf")
    rsus = load_var(flux_file, "rsus"; varname = "msuwswrf")
    rlus = load_var(flux_file, "rlus"; varname = "msuwlwrf")

    return [lhf, shf, rsus, rlus]
end

"""
    preprocess_vars(vars)

Preprocess each OutputVar in `vars` by keeping the relevant dates in
`sample_date_ranges`.
"""
function preprocess_vars(vars)
    resample_var = resampled_lonlat(CALIBRATE_CONFIG.config_file)
    vars = map(vars) do var
        var = resample_var(var)
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
    resampled_lonlat(config_file)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.
"""
function resampled_lonlat(config_file)
    config_dict = get_coupler_config_dict(CALIBRATE_CONFIG.config_file)
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
    # TODO: Generalize to 3D vars, account for stretch for 3D variables
    return var -> resampled_as(var; lon, lat)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file
    @info "Generating observations for $short_names"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    unprocessed_vars = load_vars()

    preprocessed_vars = preprocess_vars(unprocessed_vars)

    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/preprocessed_vars.jld2"),
        preprocessed_vars,
    )
    observation_vector = make_observation_vector(preprocessed_vars, sample_date_ranges)
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/obs_vec.jld2"),
        observation_vector,
    )
end
