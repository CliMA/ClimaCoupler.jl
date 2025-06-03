# TODO: Add exception for yearly vars, ensure flattening is consistent in obs map

using Distributed
import ClimaCalibrate as CAL
using ClimaCalibrate
using ClimaAnalysis
import ClimaAnalysis: SimDir, get, slice
import ClimaComms
import ClimaCoupler
using EnsembleKalmanProcesses.ParameterDistributions
using Distributions
import EnsembleKalmanProcesses as EKP
import JLD2

include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_map.jl"))
experiment_dir = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/")
model_interface = joinpath(experiment_dir, "coarse_amip", "model_interface.jl")

# Experiment Configuration
output_dir = "experiments/calibration/coarse_amip/output"
n_iterations = 9

# Prior distribution 

# constrained_gaussian(
#     "entr_param_vec",
#     [20.8, -11.6, -21.5, 16.5, -5.1, 0.024],
#     [5.0, 5.0, 5.0, 5.0, 5.0, 0.003],
#     [0.0, -30.0, -40.0, 0, -20, 0.0],
#     [40., 5, 0, 40, 15, 0.1],
# ),
# Vector of (mean, stdev) for vector parameter `entr_param_vec`
vec_distributions = [(0.0, 0.5),(0.0, 0.5),(0.0, 0.5),(0.0, 0.5),(0.0, 0.5),(0.0, 0.3), ]
vec_distributions = map(x -> Normal(x...), vec_distributions)
vec_constraints = [(0.0, 40.0), (-30.0, 5), (-40, 0), (0, 40), (-20, 15), (0.0, 0.1)]
vec_constraints = map(x -> bounded(x...), vec_constraints)
priors = [
    ParameterDistribution(VectorOfParameterized(vec_distributions), vec_constraints, "entr_param_vec"),
    constrained_gaussian("entr_mult_limiter_coeff", 1.09, 0.2, 0.0, 2.0),
    constrained_gaussian("EDMF_surface_area", 0.075, 0.02, 0.0, 0.15),
    constrained_gaussian("mixing_length_tke_surf_flux_coeff", 4, 1, 0, 8)
]
prior = combine_distributions(priors)

batch_size = 1
nt = JLD2.load_object("experiments/calibration/nt_obs.jld2")
short_names = filter(x -> x != :rsdt, keys(nt))
obs_series = create_observation_series(nt; short_names, batch_size)

eki = EKP.EnsembleKalmanProcess(obs_series, EKP.TransformUnscented(prior), verbose = true)
ensemble_size = EKP.get_N_ens(eki)

# Slurm resources for a single model run
hpc_kwargs = CAL.kwargs(time = 60 * 15, ntasks = 1, gpus_per_task = 1, cpus_per_task = 4, partition = "a3")
exeflags = "--threads=4"
eki = CAL.calibrate(CAL.GCPBackend, eki, n_iterations, prior, output_dir; model_interface, hpc_kwargs, exeflags)
