using Distributed
import ClimaCalibrate as CAL
using ClimaCalibrate
using ClimaAnalysis
import ClimaComms
import ClimaCoupler
using EnsembleKalmanProcesses.ParameterDistributions
import EnsembleKalmanProcesses as EKP
import Random
rng_seed = 1234
rng = Random.MersenneTwister(rng_seed)

include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/cld_eff_rad/observation_map.jl"))

# addprocs(CAL.SlurmManager())
# To run interactively outside of a slurm job, comment out line above and run:
# 9 works (slurm tasks) with 1 GPU and 4 cores per task
addprocs(CAL.SlurmManager(9); cpus_per_task = 4, gpus_per_task = 1, partition = "a3", time = "08:00:00")

# Make variables and the forward model available on the worker sessions
@everywhere begin
    import ClimaCoupler
    experiment_dir = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/cld_eff_rad/")
    include(joinpath(experiment_dir, "model_interface.jl"))
end

# Experiment Configuration
output_dir = "experiments/calibration/output"
n_iterations = 10
# constrained_gaussian(mean, stdev, lower_bound, upper_bound)
priors = [
    constrained_gaussian("dust_calibration_coefficient", 0, 0.2, -2, 2),
    constrained_gaussian("seasalt_calibration_coefficient", 0, 0.2, -2, 2),
    constrained_gaussian("ammonium_sulfate_calibration_coefficient", 0, 0.2, -2, 2),
    constrained_gaussian("liquid_water_specific_humidity_calibration_coefficent", 0, 0.2, -2, 2),
]
prior = combine_distributions(priors)
observation_path = joinpath(experiment_dir, "observations.jld2")
observations = JLD2.load_object(observation_path)

eki = EKP.EnsembleKalmanProcess(
    observations,
    EKP.TransformUnscented(prior; impose_prior = true),
    verbose = true,
);
# If TransformUnscented causes issues, just use TransformInversion:
# ensemble_size = 40
# eki = EKP.EnsembleKalmanProcess(
#     EKP.construct_initial_ensemble(rng, prior, ensemble_size),
#     observations,
#     EKP.TransformInversion(),
#     verbose = true,
# );
eki = CAL.calibrate(CAL.WorkerBackend, eki, n_iterations, prior, output_dir)
