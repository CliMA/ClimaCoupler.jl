# Generate and save experiment observations to disk
using ClimaAnalysis, JLD2, ClimaCoupler
include(joinpath(pkgdir(ClimaCoupler),"experiments/calibration/coarse_amip/observation_utils.jl"))

const obs_dir = "/home/ext_nefrathe_caltech_edu/calibration_obs"
const simdir = SimDir(joinpath(pkgdir(ClimaCoupler),"experiments/calibration/output/iteration_000/member_001/model_config/"))

diagnostic_var2d = get_monthly_averages(simdir, "rsut")
pressure = get_monthly_averages(simdir, "pfull")
diagnostic_var3d = get_monthly_averages(simdir, "ta")
diagnostic_var3d = ClimaAnalysis.Atmos.to_pressure_coordinates(diagnostic_var3d, pressure)

nt = get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)
JLD2.save_object("experiments/calibration/coarse_amip/nt_obs.jld2", nt)
nyears = 18
observation_vec = create_observation_vector(nt, nyears)
JLD2.save_object("experiments/calibration/coarse_amip/observations.jld2", observation_vec)
