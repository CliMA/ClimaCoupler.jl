ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
import ClimaCoupler
import ClimaCalibrate
import CUDA
import Dates: Date, Second
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "code_loading.jl"))
include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "subseasonal",
        "run_calibration.jl",
    ),
)
using Pkg
function ClimaCalibrate.forward_model(iter, member)
    Pkg.status()
    (; config_file, output_dir, sample_date_ranges, spinup, extend) = CALIBRATE_CONFIG
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    output_dir_root = output_dir

    # Update start date and length of simulation
    start_date = first(sample_date_ranges[iter + 1]) - spinup
    end_date = last(sample_date_ranges[iter + 1]) + extend
    ClimaCoupler.CalibrateTools.update_tspan!(config_dict, start_date, end_date)

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    ClimaCoupler.CalibrateTools.add_parameter_filepath!(config_dict, sampled_parameter_file)

    # Set member output directory
    member_output_dir =
        ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir

    @info "Simulation dates" start_date end_date
    setup_and_run(config_dict)
    @info "Completed member $member"
    return nothing
end
