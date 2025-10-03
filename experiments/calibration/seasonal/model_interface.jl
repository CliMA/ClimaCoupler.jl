ENV["CLIMACOMMS_DEVICE"] = "CUDA"
import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))
const config_file = joinpath(pkgdir(ClimaCoupler), "config/subseasonal_configs/wxquest_diagedmf.yml")

function ClimaCalibrate.forward_model(iter, member)
    config_dict = get_coupler_config_dict(config_file)

    output_dir_root = config_dict["coupler_output_dir"]
    eki = ClimaCalibrate.load_ekp_struct(output_dir_root, iter)
    minibatch = EKP.get_current_minibatch(eki)
    start_date = minibatch_to_start_date(minibatch)
    config_dict["start_date"] = start_date
    @info "Current minibatch: $minibatch"
    @info "Current start date: $start_date"

    config_dict["bucket_initial_condition"] = "/glade/campaign/univ/ucit0011/cchristo/initial_conditions_v_0.5/era5_bucket_processed_$(start_date)_0000.nc"

    config_dict["t_end"] = "365days"

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
    return sim
end

function minibatch_to_start_date(batch)
    start_year = minimum(batch) + 2017
    @assert start_year >= 2018
    return "$(start_year)0901"
end

import ClimaCore: Spaces
CS() = CoupledSimulation(config_file)
