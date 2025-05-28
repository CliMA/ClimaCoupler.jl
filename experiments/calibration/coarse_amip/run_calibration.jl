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
vec_distributions = [(20.8, 5.0), (-11.6, 5.0), (-21.5, 5.0), (16.5, 5.0), (-5.1, 5.0), (0.024, 0.003)]
vec_distributions = map(x -> Normal(x...), vec_distributions)
vec_constraints = [(0.0, 40.0), (-30.0, 5), (-40, 0), (0, 40), (-20, 15), (0.0, 0.1)]
vec_constraints = map(x -> bounded(x...), vec_constraints)
priors = [
    ParameterDistribution(VectorOfParameterized(vec_distributions), vec_constraints, "entr_param_vec"),
    constrained_gaussian("entr_mult_limiter_coeff", 1.09, 0.2, 0.0, 2.0),
    constrained_gaussian("EDMF_surface_area", 0.075, 0.02, 0.05, 0.1),
    constrained_gaussian("mixing_length_tke_surf_flux_coeff", 4, 1, 0, 8)
]
prior = combine_distributions(priors)
observation_path = joinpath(experiment_dir, "observations.jld2")
observation_vec = JLD2.load_object(observation_path)

# Create the EKP.ObservationSeries
batch_size = 1
num_batches = cld(length(observation_vec), batch_size)
batches = map(1:num_batches) do i
    start_idx = (i - 1) * batch_size + 1
    end_idx = min(i * batch_size, length(observation_vec))
    start_idx:end_idx
end
minibatcher = EKP.FixedMinibatcher(batches)

series_names = string.(1:length(observation_vec))
observation_series = EKP.ObservationSeries(observation_vec, minibatcher, series_names)

eki = EKP.EnsembleKalmanProcess(observation_series, EKP.TransformUnscented(prior), verbose = true)
ensemble_size = EKP.get_N_ens(eki)

# Slurm resources for a single model run
hpc_kwargs = CAL.kwargs(time = 60 * 15, ntasks = 1, gpus_per_task = 1, cpus_per_task = 4, partition = "a3")
exeflags = "--threads=4"
eki = CAL.calibrate(CAL.GCPBackend, eki, n_iterations, prior, output_dir; model_interface, hpc_kwargs, exeflags)
