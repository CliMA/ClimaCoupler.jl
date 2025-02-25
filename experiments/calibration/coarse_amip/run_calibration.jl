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

ENV["JULIA_WORKER_TIMEOUT"] = "300.0"
addprocs(CAL.SlurmManager())
# Make variables and the forward model available on the worker sessions
@everywhere begin
    import ClimaCoupler
    experiment_dir = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/")
    include(joinpath(experiment_dir, "model_interface.jl"))
end

# Experiment Configuration
output_dir = "experiments/calibration/output"
ensemble_size = 20
n_iterations = 18 # Cycle through all data
priors = [constrained_gaussian("liquid_cloud_effective_radius", 14e-6, 6e-6, 2.5e-6, 21.5e-6)]
prior = combine_distributions(priors)
observation_vec = JLD2.load_object(joinpath(experiment_dir, "observations.jld2"))

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

eki = EKP.EnsembleKalmanProcess(
    EKP.construct_initial_ensemble(prior, ensemble_size),
    observation_series,
    EKP.TransformInversion(),
)

eki = CAL.calibrate(CAL.WorkerBackend, eki, ensemble_size, n_iterations, prior, output_dir)

# Postprocessing
import Statistics: var, mean
using Test
import CairoMakie

function scatter_plot(eki::EKP.EnsembleKalmanProcess)
    f = CairoMakie.Figure(resolution = (800, 600))
    ax = CairoMakie.Axis(f[1, 1], ylabel = "Parameter Value", xlabel = "Top of atmosphere radiative SW flux")

    g = vec.(EKP.get_g(eki; return_array = true))
    params = vec.((EKP.get_ϕ(prior, eki)))

    for (gg, uu) in zip(g, params)
        CairoMakie.scatter!(ax, gg, uu)
    end

    CairoMakie.vlines!(ax, observations, linestyle = :dash)

    output = joinpath(output_dir, "scatter.png")
    CairoMakie.save(output, f)
    return output
end

function param_versus_iter_plot(eki::EKP.EnsembleKalmanProcess)
    f = CairoMakie.Figure(resolution = (800, 600))
    ax = CairoMakie.Axis(f[1, 1], ylabel = "Parameter Value", xlabel = "Iteration")
    params = EKP.get_ϕ(prior, eki)
    for (i, param) in enumerate(params)
        CairoMakie.scatter!(ax, fill(i, length(param)), vec(param))
    end

    output = joinpath(output_dir, "param_vs_iter.png")
    CairoMakie.save(output, f)
    return output
end

scatter_plot(eki)
param_versus_iter_plot(eki)

params = EKP.get_ϕ(prior, eki)
spread = map(var, params)

# Spread should be heavily decreased as particles have converged
@test last(spread) / first(spread) < 0.1
