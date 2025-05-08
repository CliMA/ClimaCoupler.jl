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
output_dir = "experiments/calibration/coarse_amip/output_simple"
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
nt = JLD2.load_object("experiments/calibration/nt_obs_3d_8_h_elem.jld2")
short_names = [:rsut, :sw_cre, :lwp]
model_error_scale = 0.10
regularization = 1e-3
obs_series = create_observation_series(nt; model_error_scale, regularization, short_names, batch_size)

# module_load_str = ClimaCalibrate.module_load_str(get_backend())
# ClimaCalibrate.slurm_model_run(0, 1, output_dir, experiment_dir, model_interface; module_load_str, exeflags)
eki = EKP.EnsembleKalmanProcess(obs_series, EKP.TransformUnscented(prior, impose_prior = true), verbose = true)
ensemble_size = EKP.get_N_ens(eki)

# Slurm resources for a single model run
hpc_kwargs = CAL.kwargs(time = 60 * 12, ntasks = 2, gpus_per_task = 1, cpus_per_task = 4, M = "nat.henrici@gmail.com", m = "ea")
exeflags = "--threads=4"
reruns = 2
eki = CAL.calibrate(CAL.DerechoBackend, eki, n_iterations, prior, output_dir; model_interface, hpc_kwargs, exeflags, reruns)
