using Distributed
import ClimaCalibrate as CAL
using ClimaCalibrate
using ClimaAnalysis
import ClimaComms
import ClimaCoupler
using EnsembleKalmanProcesses.ParameterDistributions
import EnsembleKalmanProcesses as EKP
import Random
import JLD2
rng_seed = 1234
rng = Random.MersenneTwister(rng_seed)

FULL_CALIBRATION = false

include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/cld_eff_rad/observation_map.jl"))

addprocs(CAL.SlurmManager())
# To run interactively outside of a slurm job, comment out line above and run:
# 9 works (slurm tasks) with 1 GPU and 4 cores per task
# addprocs(CAL.SlurmManager(9); cpus_per_task = 4, gpus_per_task = 1, partition = "a3", time = "08:00:00")

# Make variables and the forward model available on the worker sessions
@everywhere begin
    import ClimaCoupler
    experiment_dir = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/cld_eff_rad/")
    include(joinpath(experiment_dir, "model_interface.jl"))
end

# Experiment Configuration
output_dir = "experiments/calibration/output"
# isdir(output_dir) && error("$output_dir already exists!")
n_iterations = 10
# constrained_gaussian(mean, stdev, lower_bound, upper_bound)
# priors = [
#     # constrained_gaussian("dust_calibration_coefficient", 0, 0.2, -2, 2),
#     # constrained_gaussian("seasalt_calibration_coefficient", 0, 0.2, -2, 2),
#     # constrained_gaussian("ammonium_sulfate_calibration_coefficient", 0, 0.2, -2, 2),
#     constrained_gaussian("prescribed_cloud_droplet_number_concentration", 3e8, 5e7, 1e7, 1e9), # N₀
#     constrained_gaussian("reference_liquid_water_specific_humidity", 5e-6, 1e-6, 1e-7, 1e-5), # q₀_liq
#     constrained_gaussian("liquid_water_specific_humidity_calibration_coefficent", 0.0, 1.0, -3, 3) # α_q_liq
#     # constrained_gaussian("ice_cloud_effective_radius", 30e-6, 5e-6, 5e-6, 60e-6)
# ]

priors = [
    constrained_gaussian("ice_cloud_effective_radius", 25e-6, 1e-5, 5e-6, 9e-5), # r_eff ice
    constrained_gaussian("prescribed_cloud_droplet_number_concentration", 1e8, 5e7, 1e7, 1e9), # q₀_liq
    constrained_gaussian("reference_liquid_water_specific_humidity", 0.1, 0.05, 0.01, 0.3), # q₀_liq
    constrained_gaussian("liquid_water_specific_humidity_calibration_coefficent", 0.0, 1.0, -3.0, 3.0) # α_q_liq
]

# priors = [
# constrained_gaussian("prescribed_cloud_droplet_number_concentration", 1e8, 5e7, 1e7, 1e9), # N₀_liq
# constrained_gaussian("reference_liquid_water_specific_humidity", 0.1, 0.05, 0.01, 0.3), # q₀_liq
# constrained_gaussian("liquid_water_specific_humidity_calibration_coefficent", 0.0, 1.0, -3, 3), # α_q_liq
# constrained_gaussian("dust_calibration_coefficient", 0.0, 1.0, -3, 3), # α_dust,
# constrained_gaussian("seasalt_calibration_coefficient", 0.0, 1.0, -3, 3), # α_seasalt,
# constrained_gaussian("ammonium_sulfate_calibration_coefficient", 0.0, 1.0, -3, 3), # α_SO4,
# constrained_gaussian("reference_dust_aerosol_mass_concentration", 1e-8, 4e-9, 1e-9, 1e-6), # c₀_dust
# constrained_gaussian("reference_seasalt_aerosol_mass_concentration", 1e-8, 4e-9, 1e-9, 1e-6), # c₀_seasalt
# constrained_gaussian("reference_ammonium_sulfate_mass_concentration", 1e-8, 4e-9, 1e-9, 1e-6) # c₀_SO4
# ]

prior = combine_distributions(priors)
observation_path = joinpath(experiment_dir, "sep_obs_series.jld2")
observations = JLD2.load_object(observation_path)

eki = EKP.EnsembleKalmanProcess(
    observations,
    EKP.TransformUnscented(prior; impose_prior = true),
    verbose = true,
    scheduler = EKP.DataMisfitController(terminate_at=100),
);
ensemble_size = EKP.get_N_ens(eki)
# If TransformUnscented causes issues, just use TransformInversion:
# ensemble_size = 20
# eki = EKP.EnsembleKalmanProcess(
#     EKP.construct_initial_ensemble(rng, prior, ensemble_size),
#     observations,
#     EKP.TransformInversion(),
#     verbose = true,
#     scheduler = EKP.DataMisfitController(terminate_at=100)
# );
eki = CAL.calibrate(CAL.WorkerBackend, eki, n_iterations, prior, output_dir)
