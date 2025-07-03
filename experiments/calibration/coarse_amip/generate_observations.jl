# Generate and save experiment observations to disk
using ClimaAnalysis, JLD2, ClimaCoupler
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_utils.jl"))

const obs_dir = "/home/ext_nefrathe_caltech_edu/calibration_obs"
const simdir = SimDir(
    joinpath(
        pkgdir(ClimaCoupler),
        "output/output_0000",
    ),
)

diagnostic_var2d = get_monthly_averages(simdir, "rsut")
pressure = get_monthly_averages(simdir, "pfull")

min_z = minimum(filter(x -> x > 80, pressure.dims[altitude_name(pressure)]))
pressure = window(pressure, "z", left = min_z)
diagnostic_var3d = get_monthly_averages(simdir, "ta")
diagnostic_var3d = window(diagnostic_var3d, "z", left = min_z)
diagnostic_var3d = ClimaAnalysis.Atmos.to_pressure_coordinates(diagnostic_var3d, pressure)

nt = get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)
JLD2.save_object("experiments/calibration/nt_obs_3d_32_h_elem.jld2", nt)
