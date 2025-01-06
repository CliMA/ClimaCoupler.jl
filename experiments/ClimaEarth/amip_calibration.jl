using Distributed
import ClimaCalibrate as CAL
using ClimaCalibrate
import ClimaAnalysis: SimDir, get, slice, average_xy

include(joinpath(CAL.project_dir(), "calibration_interface.jl"))
CAL.forward_model(0,1)
