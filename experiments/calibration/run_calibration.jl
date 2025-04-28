using Distributed
import ClimaCalibrate as CAL
using ClimaCalibrate
using ClimaAnalysis
import ClimaAnalysis: SimDir, get, slice, average_xy
import ClimaComms
import EnsembleKalmanProcesses: I, ParameterDistributions.constrained_gaussian
import EnsembleKalmanProcesses as EKP
using Statistics

# Ensure ClimaComms doesn't use MPI
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
ClimaComms.@import_required_backends

single_member_dims = 1
function CAL.observation_map(iteration)
    G_ensemble = Array{Float64}(undef, single_member_dims, ensemble_size)

    for m in 1:ensemble_size
        member_path = CAL.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "model_config/output_active/clima_atmos")
        if isdir(simdir_path)
            simdir = SimDir(simdir_path)
            G_ensemble[:, m] .= process_member_data(simdir)
        else
            @info "No data found for member $m."
            G_ensemble[:, m] .= NaN
        end
    end
    return G_ensemble
end

function process_member_data(simdir::SimDir)
    output = zeros(single_member_dims)
    days = 86_400

    rsut = get(simdir; short_name = "rsut", reduction = "average", period = "30d")
    rsut_slice = slice(average_lon(average_lat(rsut)); time = 30days).data
    return rsut_slice
end

addprocs(CAL.SlurmManager())
# Make variables and the forward model available on the worker sessions
@everywhere begin
    import ClimaComms, CUDA
    ENV["CLIMACOMMS_DEVICE"] = "CUDA"
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
    import ClimaCalibrate as CAL
    import JLD2

    experiment_dir = CAL.project_dir()
    include(joinpath(experiment_dir, "model_interface.jl"))
    output_dir = joinpath(experiment_dir, "output")
    obs_path = joinpath(experiment_dir, "observations.jld2")
end

# Experiment Configuration
n_iterations = 5
noise = 0.1 * I
prior = constrained_gaussian("total_solar_irradiance", 1000, 500, 250, 2000)

# Generate observations if needed
if !isfile(obs_path)
    import JLD2
    @info "Generating observations"
    obs_output_dir = CAL.path_to_ensemble_member(output_dir, 0, 0)
    mkpath(obs_output_dir)
    touch(joinpath(obs_output_dir, "parameters.toml"))
    CAL.forward_model(0, 0)
    observations = Vector{Float64}(undef, 1)
    observations .= process_member_data(SimDir(joinpath(obs_output_dir, "model_config/output_active/clima_atmos")))
    JLD2.save_object(obs_path, observations)
end
observations = JLD2.load_object(obs_path)
@show observations
u0_mean = [mean(prior)]
uu0_cov = cov(prior)
eki = EKP.EnsembleKalmanProcess(
    EKP.ObservationSeries(EKP.Observation(observations, noise, "rsut")),
    EKP.Unscented(u0_mean, uu0_cov),
)
ensemble_size = EKP.get_N_ens(eki)

eki = CAL.calibrate(CAL.WorkerBackend, eki, n_iterations, prior, output_dir)

# Postprocessing
import EnsembleKalmanProcesses as EKP
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
@test last(spread) / first(spread) < 0.15
