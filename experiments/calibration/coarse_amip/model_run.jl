ENV["CLIMACOMMS_DEVICE"] = "CUDA"
import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))
const config_file = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "coarse_amip", "model_config.yml")

config_dict = get_coupler_config_dict(config_file)
println("tmpdir")
println(ENV["TMPDIR"])

output_dir_root = config_dict["coupler_output_dir"]
config_dict["t_end"] = "10days"

sim = try
    # Ensure that the most recent `setup_and_run` method is used, preventing
    # world age errors.
    Base.invokelatest(setup_and_run, config_dict)
catch e
    @error e
    println(catch_backtrace())
    rethrow(e)
end
