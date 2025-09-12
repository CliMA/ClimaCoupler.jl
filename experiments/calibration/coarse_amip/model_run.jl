ENV["CLIMACOMMS_DEVICE"] = "CUDA"
import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))
const config_file = joinpath(pkgdir(ClimaCoupler), "config/subseasonal_configs/wxquest_diagedmf.yml")

config_dict = get_coupler_config_dict(config_file)
println("tmpdir")
println(ENV["TMPDIR"])
config_dict["t_end"] = "366days"
@show config_dict["coupler_output_dir"]
Base.invokelatest(setup_and_run, config_dict)

