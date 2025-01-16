using Distributed
import ClimaCalibrate as CAL
using ClimaCalibrate
using ClimaAnalysis
import ClimaAnalysis: SimDir, get, slice, average_xy
import ClimaComms
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
ClimaComms.@import_required_backends

function CAL.observation_map(iteration)
    single_member_dims = (2,)
    G_ensemble = Array{Float64}(undef, single_member_dims..., ensemble_size)

    for m in 1:ensemble_size
        member_path = CAL.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "output_active")
        if isdir(simdir_path)
            simdir = SimDir(simdir_path)
            G_ensemble[:, m] .= process_member_data(simdir)
        else
            G_ensemble[:, m] .= NaN
        end
    end
    return G_ensemble
end

function process_member_data(simdir::SimDir)
    isempty(simdir) && return NaN
    rsut =
    try
        get(simdir; short_name = "rsut", reduction = "average", period = "30d")
    catch e
        @error e
        return NaN
    end
    return slice(average_lon(average_lat(rsut)); time = 30).data
end

# addprocs(CAL.SlurmManager(5))

@everywhere begin
    import ClimaComms, CUDA
    ENV["CLIMACOMMS_DEVICE"] = "CUDA"
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
    import ClimaCalibrate as CAL
    import ClimaAtmos as CA
    import JLD2
    import EnsembleKalmanProcesses:
        I, ParameterDistributions.constrained_gaussian, ParameterDistributions.combine_distributions

    experiment_dir = CAL.project_dir()
    include(joinpath(experiment_dir, "calibration_interface.jl"))
    output_dir = joinpath(experiment_dir, "output")

    # Experiment Configuration
    ensemble_size = 50
    n_iterations = 5
    astronomical_unit = 149_597_870_000
    noise = 0.1 * I
    priors = [constrained_gaussian("astronomical_unit", 1.5e11, 1e11, 2e5, Inf),
[entr_inv_tau] #0.001 - 0.004

            constrained_gaussian("specific_humidity_precipitation_threshold", 5.0e-6, 5.0e-6,  1e-9, Inf),
            constrained_gaussian("precipitation_timescale", 500, 300, 240, 1000),
            # constrained_gaussian("supersaturation_precipitation_threshold", 0.04, 0.03, 0, 0.1),
        ]
    prior = combine_distributions(priors)
    obs_path = joinpath(experiment_dir, "observations.jld2")
end

# Generate observations if needed
if !isfile(obs_path)
    import JLD2
    @info "Generating observations"
    obs_output_dir = CAL.path_to_ensemble_member(output_dir, 0, 0)
    mkpath(obs_output_dir)
    touch(joinpath(obs_output_dir, "parameters.toml"))
    @everywhere begin
        CAL.forward_model(0, 0)
    end
    observations = Vector{Float64}(undef, 1)
    observations .= process_member_data(SimDir(joinpath(obs_output_dir, "output_active")))
    JLD2.save_object(obs_path, observations)
end

# Initialize experiment data
@everywhere observations = JLD2.load_object(obs_path)

eki = CAL.calibrate(
    CAL.WorkerBackend,
    ensemble_size,
    n_iterations,
    observations,
    noise,
    prior,
    output_dir,
)

# TODO: Enable `calibrate` to checkpoint, rerunning from midway through calibration
# Postprocessing
import EnsembleKalmanProcesses as EKP
import Statistics: var, mean
using Test
import CairoMakie

function scatter_plot(eki::EKP.EnsembleKalmanProcess)
    f = CairoMakie.Figure(resolution = (800, 600))
    ax = CairoMakie.Axis(
        f[1, 1],
        ylabel = "Parameter Value",
        xlabel = "Top of atmosphere radiative SW flux",
    )

    g = vec.(EKP.get_g(eki; return_array = true))
    params = vec.((EKP.get_ϕ(prior, eki)))

    for (gg, uu) in zip(g, params)
        CairoMakie.scatter!(ax, gg, uu)
    end

    CairoMakie.hlines!(ax, [astronomical_unit], linestyle = :dash)
    CairoMakie.vlines!(ax, observations, linestyle = :dash)

    output = joinpath(output_dir, "scatter.png")
    CairoMakie.save(output, f)
    return output
end

function param_versus_iter_plot(eki::EKP.EnsembleKalmanProcess)
    f = CairoMakie.Figure(resolution = (800, 600))
    ax = CairoMakie.Axis(
        f[1, 1],
        ylabel = "Parameter Value",
        xlabel = "Iteration",
    )
    params = EKP.get_ϕ(prior, eki)
    for (i, param) in enumerate(params)
        CairoMakie.scatter!(ax, fill(i, length(param)), vec(param))
    end

    CairoMakie.hlines!(ax, [astronomical_unit]; color = :red, linestyle = :dash)

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
# Parameter should be close to true value
@test mean(last(params)) ≈ astronomical_unit rtol = 0.02
