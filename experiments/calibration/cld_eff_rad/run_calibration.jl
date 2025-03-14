using Distributed
import ClimaCalibrate as CAL
using ClimaCalibrate
using ClimaAnalysis
import ClimaAnalysis: SimDir, get, slice, average_xy
import ClimaComms
import ClimaCoupler
using EnsembleKalmanProcesses.ParameterDistributions
import EnsembleKalmanProcesses as EKP
import Random
rng_seed = 1234
rng = Random.MersenneTwister(rng_seed)

include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/cld_eff_rad/observation_map.jl"))

addprocs(CAL.SlurmManager())

# Make variables and the forward model available on the worker sessions
@everywhere begin
    import ClimaCoupler
    experiment_dir = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/cld_eff_rad/")
    include(joinpath(experiment_dir, "model_interface.jl"))
end

# Experiment Configuration
output_dir = "experiments/calibration/output"
ensemble_size = 10
n_iterations = 10
priors = [
    constrained_gaussian("prescribed_cloud_droplet_number_concentration", 1.5e8, 9e7, 1e7, 1e9),
]
prior = combine_distributions(priors)
obs_path = joinpath(experiment_dir, "observations.jld2")
observations = JLD2.load_object(obs_path)

eki = EKP.EnsembleKalmanProcess(
    EKP.construct_initial_ensemble(rng, prior, ensemble_size),
    observations,
    EKP.TransformInversion(),
    verbose=true
)

eki = CAL.calibrate(CAL.WorkerBackend, eki, ensemble_size, n_iterations, prior, output_dir)
