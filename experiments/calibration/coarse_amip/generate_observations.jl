# Generate and save experiment observations to disk
using ClimaAnalysis, JLD2
include("experiments/calibration/coarse_amip/observation_utils.jl")

const obs_dir = "/home/ext_nefrathe_caltech_edu/calibration_obs"
const diagnostic_dir = "experiments/calibration/output/old/iteration_000/member_001/model_config/output_0000/clima_atmos/"

diagnostic_var2d = OutputVar(joinpath(diagnostic_dir, "rsdt_1M_average.nc"));
pressure = OutputVar(joinpath(diagnostic_dir, "pfull_1M_average.nc"));
diagnostic_var3d = OutputVar(joinpath(diagnostic_dir, "ta_1M_average.nc"));
diagnostic_var3d = ClimaAnalysis.Atmos.to_pressure_coordinates(diagnostic_var3d, pressure)

nt = get_all_output_vars(obs_dir, diagnostic_var2d, diagnostic_var3d)
nyears = 18
observation_vec = create_observation_vector(nt, nyears)
JLD2.save_object("experiments/calibration/coarse_amip/observations.jld2", observation_vec)
