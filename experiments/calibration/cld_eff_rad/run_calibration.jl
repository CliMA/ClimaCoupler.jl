using Distributed
import ClimaCalibrate as CAL
using ClimaCalibrate
using ClimaAnalysis
import ClimaAnalysis: SimDir, get, slice, average_xy
import ClimaComms
import ClimaCoupler
using EnsembleKalmanProcesses.ParameterDistributions
import EnsembleKalmanProcesses as EKP

include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_map.jl"))

addprocs(CAL.SlurmManager())

# Make variables and the forward model available on the worker sessions
@everywhere begin
    import ClimaCoupler
    experiment_dir = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/cld_eff_rad/")
    include(joinpath(experiment_dir, "model_interface.jl"))
end

# Experiment Configuration
output_dir = "experiments/calibration/output"
ensemble_size = 30
n_iterations = 9
priors = [
    constrained_gaussian("liquid_cloud_effective_radius", 14e-6, 6e-6, 2.5e-6, 21.5e-6),
    constrained_gaussian("ice_cloud_effective_radius", 25e-6, 6e-6, 2.5e-6, 33e-6),
]
prior = combine_distributions(priors)

# TODO: Add observation

eki = EKP.EnsembleKalmanProcess(
    EKP.construct_initial_ensemble(prior, ensemble_size),
    observation_series,
    EKP.TransformInversion(),
    verbose=true
)

eki = CAL.calibrate(CAL.WorkerBackend, eki, ensemble_size, n_iterations, prior, output_dir)
