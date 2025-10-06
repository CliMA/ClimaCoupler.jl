import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2
import ClimaDiagnostics
import ClimaCore
# Access CalibrateConfig
include(joinpath(@__DIR__, "run_calibration.jl"))

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/data_sources.jl"))

function make_observation_vector(
    short_names,
    sample_date_ranges,
    obsdir,
    config_file,
)
    vars = load_vars(obsdir, short_names)
    vars = preprocess_vars(vars, sample_date_ranges, config_file)
    return make_observation_vector(vars, sample_date_ranges)
end

set_short_name!(var, short_name) = (var.attributes["short_name"] = short_name)

"""
    load_vars(obsdir, short_names)
Load NetCDF files belonging to `short_names` in `obsdir` as `OutputVar`s.
"""
function load_vars(obsdir, short_names)
    rad_and_pr_obs_dict = get_obs_var_dict()
    rsut = rad_and_pr_obs_dict["rsut"](start_date)
    rsutcs = rad_and_pr_obs_dict["rsutcs"](start_date)
    sw_cre = rsutcs - rsut
    set_short_name!(sw_cre, "sw_cre")
    return [sw_cre]
end

"""
    preprocess_vars(vars, sample_date_ranges, config_file)
Preprocess each OutputVar in `vars` by keeping the relevant dates in
`sample_date_ranges`.
"""
function preprocess_vars(vars, sample_date_ranges, config_file)
    out_var = OutputVar("/glade/derecho/scratch/nefrathe/tmp/output_seasonal/iteration_000/member_001/amip_config/output_0000/clima_atmos/rsut_1M_average.nc")
    resample_var(x) = ClimaAnalysis.resampled_as(x, out_var ; dim_names = ["lon", "lat"])
    vars = resample_var.(vars)
    vars = map(vars) do var
        set_units(var, "W m^-2")
    end
    return vars
end

var_units = Dict(
    "pr" => "kg m^-2 s^-1",
    "mslp" => "Pa",
    "tas" => "K",
    "sw_cre" => "W m^-2"
)


function make_observation_vector(vars, covariance_sample_ranges)
    obs_vec = map(CALIBRATE_CONFIG.sample_date_ranges) do sample_date_range
        covar_estimator = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(
            covariance_sample_ranges;
            model_error_scale = 0.05,
            regularization = 1e-5,
        )
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
    ENV["CLIMACOMMS_DEVICE"] = "CPU"

    cs = CoupledSimulation(config_file)
    ENV["CLIMACOMMS_DEVICE"] = "CUDA"
    center_space = cs.model_sims.atmos_sim.domain.center_space
    (lon_nlevels, lat_nlevels, z_nlevels) =
        ClimaDiagnostics.Writers.default_num_points(center_space)
    longitudes = range(-180, 180, lon_nlevels)
    latitudes = range(-90, 90, lat_nlevels)
    stretch = center_space.grid.vertical_grid.topology.mesh.stretch
    # TODO: Account for stretch for 3D variables and interpolate to pressure?
    dz_bottom = center_space.grid.vertical_grid.topology.mesh.faces[2].z
    z_levels = range(dz_bottom, ClimaCore.Spaces.z_max(center_space), z_nlevels)
    return var -> resampled_as(var; lon = longitudes, lat = latitudes)
end


if abspath(PROGRAM_FILE) == @__FILE__
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file
    @info "Generating observations for $short_names"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    # TODO: Use resample func
    start_date = DateTime(2009,9,1)

    unprocessed_vars = load_vars(nothing, short_names)
    preprocessed_vars = preprocess_vars(unprocessed_vars, sample_date_ranges, config_file)
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler),"experiments/calibration/era5_preprocessed_vars.jld2"),
        preprocessed_vars,
    )
    covariance_sample_ranges = [(DateTime(yr, 1, 1), DateTime(yr, 3, 1)) for yr in 2001:2019]
    observation_vector = make_observation_vector(preprocessed_vars, covariance_sample_ranges)
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler),"experiments/calibration/weatherquest_obs_vec.jld2"),
        observation_vector,
    )
end
