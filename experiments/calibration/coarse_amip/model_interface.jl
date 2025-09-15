ENV["CLIMACOMMS_DEVICE"] = "CUDA"
import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))
const config_file = joinpath(pkgdir(ClimaCoupler), "config/subseasonal_configs/wxquest_diagedmf.yml")

function ClimaCalibrate.forward_model(iter, member)

    config_dict = get_coupler_config_dict(config_file)
    println("tmpdir")
    println(ENV["TMPDIR"])

    output_dir_root = config_dict["coupler_output_dir"]
    eki = ClimaCalibrate.load_ekp_struct(output_dir_root, iter)
    minibatch = EKP.get_current_minibatch(eki)
    
    config_dict["start_date"] = minibatch_to_start_date(minibatch)
    @info "Current minibatch: $minibatch"
    @info "Current start date: $(minibatch_to_start_date(minibatch))"

    config_dict["bucket_initial_condition"] = "/glade/campaign/univ/ucit0011/cchristo/wxquest_ics/era5_bucket_processed_$(start_date)_0000.nc"

    config_dict["t_end"] = "39days"

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    config_dict["coupler_toml"] = [sampled_parameter_file]
    # Set member output directory
    member_output_dir = ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir

    sim = try
        # Ensure that the most recent `setup_and_run` method is used, preventing
        # world age errors.
        Base.invokelatest(setup_and_run, config_dict)
    catch e
        @error e
        println(catch_backtrace())
        rethrow(e)
    end

    @info "Completed member $member"
    @info "Removing checkpoints"
    # rm(joinpath(member_output_dir, "model_config", "checkpoints"), recursive = true)
    return sim
end

function minibatch_to_start_date(batch)
    start_year = minimum(batch) + 2017
    @assert start_year >= 2018
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
