import ClimaAnalysis
import CairoMakie
import EnsembleKalmanProcesses as EKP

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
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
end

function plot_g_ens_mean(output_dir, ekp)

    var = EKP.get_observation_series(ekp) |> EKP.get_metadata
    size_tuple = size(var.data)
    G_ens_mean = reshape(get_g_mean_final(ekp), size_tuple)
    var = remake(var, data = G_ens_mean)

    for (i, t) in enumerate(dates(var))
        sliced_var = slice(var, time = t)

        # TODO: Make it easier to adjust the colorrange for more variables
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[1, i], sliced_var, :plot => Dict(
            :title => "$(short_name(var)) g ensemble mean for\n$t"
            :colorrange => (-200.0, 50.0),
        ))
    end
end

function plot_pointwise_loss(output_dir, iter, ekp)
end

function plot_obs(output_dir, iter, ekp)
    iter > 0 && return

    var = EKP.get_observation_series(ekp) |> EKP.get_metadata
    size_tuple = size(var.data)
    data = reshape(get_obs(ekp), size_tuple)
    var = remake(var, data = data)

    # TODO: Make this into a square instead of a long strip
    fig = Figure(length(times(var)), size = (length(times(var)) * 500, 500))

    for (i, t) in enumerate(dates(var))
        sliced_var = slice(var, time = t)

        # TODO: Make it easier to adjust the colorrange for more variables
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[1, i], sliced_var, :plot => Dict(
            :title => "$(short_name(var)) observational data for\n$t"
            :colorrange => (-200.0, 50.0),
        ))
    end
end

function plot_obs_noise_cov(output_dir, iter, ekp)
    iter > 0 && return

    var = EKP.get_observation_series(ekp) |> EKP.get_metadata
    size_tuple = size(var.data)
    covariance = reshape(LinearAlgebra.diag(EKP.get_obs_noise_cov(ekp,build=false)[1]), size_tuple)
    var = remake(var, data = covariance)

    # TODO: Make this into a square instead of a long strip
    fig = Figure(size = (length(times(var)) * 800, 600))

    for (i, t) in enumerate(dates(var))
        sliced_var = slice(var, time = t)
        max_val = round(maximum(sliced_var.data), sigdigits = 3)
        min_val = round(minimum(sliced_var.data), sigdigits = 3)
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[1, i], sliced_var, more_kwargs = Dict(
            :axis => Dict(:title => "Variance of $(short_name(var)) for $t\n[minimum: $min_val and maximum: $max_val]"),
            :plot => Dict(
                :colorrange => (0.0, 1000.0),
                :colormap => CairoMakie.cgrad(:viridis, 21, categorical = true),
            ),
        ))
    end
    save(joinpath(output_dir, "variance_heatmap.png"), fig)
end
