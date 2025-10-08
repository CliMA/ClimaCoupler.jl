ENV["CLIMACOMMS_DEVICE"] = "CUDA"
import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))
include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "run_calibration.jl"))

function ClimaCalibrate.forward_model(iter, member)

    config_dict = get_coupler_config_dict(CALIBRATE_CONFIG.config_file)
    output_dir_root = CALIBRATE_CONFIG.output_dir
    eki = ClimaCalibrate.load_ekp_struct(output_dir_root, iter)
    minibatch = EKP.get_current_minibatch(eki)
    start_date = first(CALIBRATE_CONFIG.sample_date_ranges[iter+1]) - CALIBRATE_CONFIG.spinup
    start_date_str = replace(string(Date(start_date)), "-" => "")
    end_date = last(CALIBRATE_CONFIG.sample_date_ranges[iter+1]) + CALIBRATE_CONFIG.extend
    sim_length = Second(end_date - start_date)

    config_dict["start_date"] = start_date_str
    config_dict["bucket_initial_condition"] = "/glade/campaign/univ/ucit0011/cchristo/initial_conditions_v_0.5/era5_bucket_processed_$(start_date_str)_0000.nc"
    config_dict["t_end"] = "$(sim_length)secs"

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    config_dict["coupler_toml"] = [sampled_parameter_file]
    # Set member output directory
    member_output_dir = ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir

    @info "Current minibatch: $minibatch"
    @info "Simulation dates" start_date end_date

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
