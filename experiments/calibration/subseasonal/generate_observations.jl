import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import JLD2
import ClimaDiagnostics
import ClimaCore
# Access CalibrateConfig
include(joinpath(@__DIR__, "run_calibration.jl"))

"""
    regrid_vars(vars)

Regrid each OutputVar in `vars` to match the grid of the diagnostics output from
ClimaAtmos.
"""
function regrid_vars(vars)
    resample_var = resampled_lonlat(CALIBRATE_CONFIG.config_file)
    vars = map(vars) do var
        var = resample_var(var)
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
    # TODO: Generalize to 3D vars, account for stretch for 3D variables
    return var -> ClimaAnalysis.resampled_as(var; lon, lat)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
    (; sample_date_ranges, short_names, config_file) = CALIBRATE_CONFIG
    @info "Generating observations for $short_names"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    data_loader = CalibrationTools.ERA5DataLoader()
    varnames = ["hfls", "hfss", "rsus", "rlus"]
    unprocessed_vars = get.(Ref(data_loader), varnames)
    preprocessed_vars = regrid_vars(unprocessed_vars)

    JLD2.save_object(
        joinpath(
            pkgdir(ClimaCoupler),
            "experiments/calibration/subseasonal/preprocessed_vars.jld2",
        ),
        preprocessed_vars,
    )
    observation_vector = make_observation_vector(preprocessed_vars, sample_date_ranges)
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/obs_vec.jld2"),
        observation_vector,
    )
end
