# TODO: Make exceptional variables easy to process: non-seasonal, missing data
# ERROR: 2011-09-01T00:00:00 is not a date in the OutputVar with the short name cl
# ERROR: 2017-12-01T00:00:00 is not a date in the OutputVar with the short name cl

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
output_dir = joinpath(ENV["TMPDIR"],"output_quick")
mkpath(output_dir)
n_iterations = 9

# Prior distribution 
priors = [
    constrained_gaussian("mixing_length_diss_coeff", 0.22, 0.07, 0, 1),
    # constrained_gaussian("entr_pi_const", 0.7319794747190422, 0.15, 0, 10),
    constrained_gaussian("precipitation_timescale", 919.3827604731249, 150.0, 0, Inf),
    constrained_gaussian("EDMF_surface_area", 0.10928882001604676, 0.03, 0, Inf),
]
prior = combine_distributions(priors)

batch_size = 1
nt = JLD2.load_object("/glade/u/home/nefrathe/clima/ClimaCoupler.jl/experiments/calibration/nt_obs_3d_8_h_elem.jld2")
short_names = [:rsut, :sw_cre, :lwp]
model_error_scale = 0.10
regularization = 1e-3
# obs_series = create_observation_series(nt; model_error_scale, regularization, short_names, batch_size)

eki = EKP.EnsembleKalmanProcess(EKP.Observation([1], EKP.I, "placeholder"), EKP.TransformUnscented(prior, impose_prior = true), verbose = true)
ensemble_size = EKP.get_N_ens(eki)
# Slurm resources for a single model run
hpc_kwargs = CAL.kwargs(time = 660, ntasks = 8, gpus_per_task = 1, cpus_per_task = 4, job_priority = "premium", q = "main")
exeflags = "--threads=4"
eki = CAL.calibrate(CAL.DerechoBackend, eki, n_iterations, prior, output_dir; model_interface, hpc_kwargs, exeflags)
