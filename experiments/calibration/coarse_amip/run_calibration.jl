# TODO: Add exception for yearly vars

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
output_dir = "experiments/calibration/coarse_amip/output_32_h_elem"
n_iterations = 9

# Prior distribution 
priors = [
    constrained_gaussian("pi_groups_coeff", 1.0, 0.3, 0, Inf),
    constrained_gaussian("entr_pi_const", 0.7319794747190422, 0.15, 0, 10),
    constrained_gaussian("mixing_length_tke_surf_flux_coeff", 9.233650392728526, 1.5, 0, Inf),
    constrained_gaussian("precipitation_timescale", 919.3827604731249, 150.0, 0, Inf),
    constrained_gaussian("EDMF_surface_area", 0.10928882001604676, 0.03, 0, Inf)
]
prior = combine_distributions(priors)

batch_size = 1
nt = JLD2.load_object("experiments/calibration/nt_obs_32_h_elem.jld2")
short_names = filter(x -> !(x ∈ [:rsdt, :net_rad]), keys(nt))
model_error_scale = 0.10
regularization = 1e-3
obs_series = create_observation_series(nt;model_error_scale, regularization, short_names, batch_size)

eki = EKP.EnsembleKalmanProcess(obs_series, EKP.TransformUnscented(prior, impose_prior = true), verbose = true)
ensemble_size = EKP.get_N_ens(eki)

# Slurm resources for a single model run
hpc_kwargs = CAL.kwargs(time = 60 * 36, ntasks = 4, gpus_per_task = 1, cpus_per_task = 4, partition = "a3")
exeflags = "--threads=4"
eki = CAL.calibrate(CAL.GCPBackend, eki, n_iterations, prior, output_dir; model_interface, hpc_kwargs, exeflags)
