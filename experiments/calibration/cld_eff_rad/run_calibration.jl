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

ensemble_size = 40
addprocs(CAL.SlurmManager())
# For running interactively:
# addprocs(CAL.SlurmManager(ensemble_size); cpus_per_task = 4, gpus_per_task = 1, partition = "a3", time = "08:00:00")

# Make variables and the forward model available on the worker sessions
@everywhere begin
    import ClimaCoupler
    experiment_dir = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/cld_eff_rad/")
    include(joinpath(experiment_dir, "model_interface.jl"))
end

# Experiment Configuration
output_dir = "experiments/calibration/output"
n_iterations = 10
priors = [
    constrained_gaussian("dust_calibration_coefficient", 0, 0.2, -2, 2),
    constrained_gaussian("seasalt_calibration_coefficient", 0, 0.2, -2, 2),
    constrained_gaussian("ammonium_sulfate_calibration_coefficient", 0, 0.2, -2, 2),
    constrained_gaussian("liquid_water_specific_humidity_calibration_coefficent", 0, 0.2, -2, 2),
    constrained_gaussian("reference_dust_aerosol_mass_concentration", 1e-8, 0.2, -2, 2),
    constrained_gaussian("reference_seasalt_aerosol_mass_concentration", 1e-8, 0.2, -2, 2),
    constrained_gaussian("reference_ammonium_sulfate_mass_concentration", 1e-8, 0.2, -2, 2),
]
prior = combine_distributions(priors)
observation_path = joinpath(experiment_dir, "observations.jld2")
observations = JLD2.load_object(observation_path)

eki = EKP.EnsembleKalmanProcess(
    EKP.construct_initial_ensemble(rng, prior, ensemble_size),
    observations,
    EKP.TransformInversion(),
    verbose = true,
);
# using ClimaAtmos#main
eki = CAL.calibrate(CAL.WorkerBackend, eki, n_iterations, prior, output_dir)
