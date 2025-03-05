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
    experiment_dir = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/")
    include(joinpath(experiment_dir, "model_interface.jl"))
end

# Experiment Configuration
output_dir = "experiments/calibration/output"
ensemble_size = 30
n_iterations = 9
priors = [constrained_gaussian("liquid_cloud_effective_radius", 14e-6, 6e-6, 2.5e-6, 21.5e-6),
constrained_gaussian("ice_cloud_effective_radius", 25e-6, 6e-6, 2.5e-6, 33e-6),
constrained_gaussian("precipitation_timescale", 600, 400, 0, 1200)
]
prior = combine_distributions(priors)
observation_path = joinpath(experiment_dir, "observations.jld2")
observation_vec = JLD2.load_object(observation_path)

# Create the EKP.ObservationSeries
batch_size = 2
num_batches = cld(length(observation_vec), batch_size)
batches = map(1:num_batches) do i
    start_idx = (i - 1) * batch_size + 1
    end_idx = min(i * batch_size, length(observation_vec))
    start_idx:end_idx
end
minibatcher = EKP.FixedMinibatcher(batches)

series_names = string.(1:length(observation_vec))
observation_series = EKP.ObservationSeries(observation_vec, minibatcher, series_names)

eki = EKP.EnsembleKalmanProcess(
    EKP.construct_initial_ensemble(prior, ensemble_size),
    observation_series,
    EKP.TransformInversion(),
)

eki = CAL.calibrate(CAL.WorkerBackend, eki, ensemble_size, n_iterations, prior, output_dir)
