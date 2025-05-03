import ClimaAnalysis
import CairoMakie
import EnsembleKalmanProcesses as EKP
import GeoMakie

"""
    plot_constrained_params_and_errors(output_dir, ekp, prior)

Save a figure in `output_dir` of the contrained parameters and the error over
the number of iterations.
"""
function plot_constrained_params_and_errors(output_dir, ekp, priors)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500));
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 1], ekp)
    # EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 2], ekp, error_metric = "loss")
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
end

"""
    plot_g_ens_mean(output_dir, ekp)

Plot the mean forward model evaluation from `ekp`.
"""
function plot_g_ens_mean(output_dir, ekp)
    var = EKP.get_observation_series(ekp) |> EKP.get_metadata
    size_tuple = size(var.data)
    G_ens_mean = reshape(get_g_mean_final(ekp), size_tuple)
    var = ClimaAnalysis.remake(var, data = G_ens_mean)

    fig = CairoMakie.Figure(size = (length(ClimaAnalysis.times(var)) * 500, 500))

    for (i, t) in enumerate(dates(var))
        sliced_var = ClimaAnalysis.slice(var, time = t)

        # TODO: Make it easier to adjust the colorrange for more variables
        # TODO: This is wrong and need to be fixed
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[1, i], sliced_var, :plot => Dict(
            :title => "$(short_name(var)) g ensemble mean for\n$t"
        ))
    end
    CairoMakie.save(fig, joinpath(output_dir, "g_ens_mean_iter$iter.png"))
    return nothing
end

function plot_pointwise_loss(output_dir, iter, ekp)
end

"""
    plot_obs(output_dir, iter, ekp)

Plot the observational data from `ekp`.
"""
function plot_obs(output_dir, obs_series)
    var = EKP.get_metadata(obs_series)
    size_tuple = size(var.data)
    data = reshape(EKP.get_obs(obs_series), size_tuple)
    var = ClimaAnalysis.remake(var, data = data)
    # TODO: Make this into a square instead of a long strip
    fig = CairoMakie.Figure(size = (length(ClimaAnalysis.times(var)) * 500, 500))

    for (i, t) in enumerate(ClimaAnalysis.dates(var))
        sliced_var = ClimaAnalysis.slice(var, time = t)

        # TODO: Make it easier to adjust the colorrange for more variables
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[1, i], sliced_var, more_kwargs = Dict(:plot => Dict(
            :colorrange => (-200.0, 50.0),
            :colormap => CairoMakie.cgrad(:viridis, 11, categorical = true)
        ), :axis => Dict(:title => "$(short_name(var)) observational data for\n$t",)))
    end
    CairoMakie.save(joinpath(output_dir, "obs_from_ekp.png"), fig)
    return nothing
end

"""
    plot_obs_noise_cov(output_dir, ekp)

Plot the observational noise covariance from `ekp`.

This function assumes the covariance matrix is diagonal.
"""
function plot_obs_noise_cov(output_dir, obs_series)
    # var = EKP.get_observation_series(ekp) |> EKP.get_metadata
    var = EKP.get_metadata(obs_series)
    size_tuple = size(var.data)
    covariance = reshape(LinearAlgebra.diag(EKP.get_obs_noise_cov(obs_series,build=false)[1]), size_tuple)
    var = ClimaAnalysis.remake(var, data = covariance)

    # TODO: Make this into a square instead of a long strip
    fig = CairoMakie.Figure(size = (length(ClimaAnalysis.times(var)) * 800, 600))

    for (i, t) in enumerate(ClimaAnalysis.dates(var))
        sliced_var = ClimaAnalysis.slice(var, time = t)
        max_val = round(maximum(sliced_var.data), sigdigits = 3)
        min_val = round(minimum(sliced_var.data), sigdigits = 3)
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[1, i], sliced_var, more_kwargs = Dict(
            :axis => Dict(:title => "Variance of $(short_name(var)) for $t\n[minimum: $min_val and maximum: $max_val]"),
            :plot => Dict(
                :colorrange => (0.0, 500.0),
                :colormap => CairoMakie.cgrad(:viridis, 11, categorical = true),
            ),
        ))
    end
    CairoMakie.save(joinpath(output_dir, "variance_heatmap.png"), fig)
    return nothing
end
