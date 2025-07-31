import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
import JLD2
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))
const config_file = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "coarse_amip", "model_config.yml")

function ClimaCalibrate.forward_model(iter, member)

    config_dict = get_coupler_config_dict(config_file)

    output_dir_root = config_dict["coupler_output_dir"]
    eki = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir_root, iter))
    minibatch = EKP.get_current_minibatch(eki)
    config_dict["start_date"] = minibatch_to_start_date(minibatch)
    @info "Current minibatch: $minibatch"
    @info "Current start date: $(minibatch_to_start_date(minibatch))"
    spinup_days = 92
    nyears = length(minibatch)
    t_end_days = spinup_days + 365 * nyears
    config_dict["t_end"] = "$(t_end_days)days"

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    config_dict["coupler_toml"] = [sampled_parameter_file]
    # Set member output directory
    member_output_dir = ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir

    checkpoint_dir = joinpath(member_output_dir, "model_config", "checkpoints")
    # Use the model's own previous checkpoints if they exist
    if isdir(checkpoint_dir) && !isempty(readdir(checkpoint_dir))
        config_dict["restart_t"] = nothing
        config_dict["restart_dir"] = nothing
    end

    sim = try
        setup_and_run(config_dict)
    catch e
        @error e
        println(catch_backtrace())
        rethrow(e)
    end

    @info "Completed member $member"
    @info "Removing checkpoints"
    rm(joinpath(member_output_dir, "model_config", "checkpoints"), recursive = true)
    return sim
end

function minibatch_to_start_date(batch)
    start_year = minimum(batch) + 2005
    @assert start_year >= 2006
    return "$(start_year)0901"
end

import ClimaCore: Spaces
CS() = CoupledSimulation(config_file)

function get_resample_func()
    cs = CS()
    center_space = cs.model_sims.atmos_sim.domain.center_space
    (lon_nlevels, lat_nlevels, z_nlevels) = ClimaDiagnostics.Writers.default_num_points(center_space)
    longitudes = range(-180, 180, lon_nlevels)
    latitudes = range(-90, 90, lat_nlevels)
    stretch = center_space.grid.vertical_grid.topology.mesh.stretch
    # TODO: account for stretch for 3D variables and interpolate to pressure?
    z_levels = range(dz_bottom, Spaces.z_max(center_space), z_nlevels)
    return var -> resampled_to(var; lon = longitudes, lat = latitudes)
end
