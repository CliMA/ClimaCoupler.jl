bias_plot_extrema = Dict(
    "tas" => (-6, 6),
    "tas - ta" => (-6, 6),
    "hfls" => (-50, 50),
    "hfss" => (-25, 25),
    "rsus" => (-50, 50),
    "rlus" => (-50, 50),
    "mslp" => (-1000, 1000),
    "pr" => (-1e-4, 1e-4),
    "ta_850hPa" => (-2, 2),
    "ta_500hPa" => (-2, 2),
    "ta_200hPa" => (-2, 2),
    "hur_850hPa" => (-2, 2),
    "hur_500hPa" => (-2, 2),
    "hur_200hPa" => (-2, 2),
)

"""
    plot_bias_weekly(ekp, simdir, iteration; output_dir)

Plot bias maps comparing simulation output to ERA5 observations for all variables in
`CALIBRATE_CONFIG.short_names`. ERA5 vars are reconstructed from the EKP observation
object and denormalized to physical units when normalization is enabled.
"""
function plot_bias_weekly(ekp, simdir, iteration; output_dir = simdir.simulation_path)
    (; short_names, sample_date_ranges) = CALIBRATE_CONFIG
    sample_date_range = sample_date_ranges[iteration]
    calib_start, _ = sample_date_range

    # Reconstruct ERA5 OutputVars from the EKP observation object
    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs = ClimaCalibrate.get_observations_for_nth_iteration(obs_series, iteration)

    era5_vars =
        mapreduce(ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, minibatch_obs)

    sim_vars = load_and_preprocess_vars(simdir, short_names)

    # Match sim_vars with era5_vars by short_name
    var_pairs = []
    for sim_var in sim_vars
        sn = ClimaAnalysis.short_name(sim_var)
        era5_idx = findfirst(v -> ClimaAnalysis.short_name(v) == sn, era5_vars)
        if !isnothing(era5_idx)
            push!(var_pairs, (sim_var, era5_vars[era5_idx]))
        else
            @warn "No ERA5 data found for $sn — skipping bias plot"
        end
    end

    if isempty(var_pairs)
        @warn "No matching variable pairs found for bias plotting"
        return nothing
    end

    fig = GeoMakie.Figure(size = (2000, 500 * length(var_pairs)))
    for (i, (sim_var, era5_var)) in enumerate(var_pairs)
        sn = ClimaAnalysis.short_name(sim_var)
        sim_var_t = ClimaAnalysis.select(
            sim_var;
            by = ClimaAnalysis.MatchValue(),
            time = calib_start,
        )
        era5_var_t = ClimaAnalysis.select(
            era5_var;
            by = ClimaAnalysis.MatchValue(),
            time = calib_start,
        )
        cmap_extrema = get(bias_plot_extrema, sn, extrema(sim_var_t.data))
        try
            if ClimaAnalysis.has_pressure(sim_var_t)
                for (j, pressure) in enumerate(ClimaAnalysis.pressures(sim_var_t))
                    sim_var_t_p = ClimaAnalysis.select(
                        sim_var_t;
                        by = ClimaAnalysis.MatchValue(),
                        pressure,
                    )
                    era5_var_t_p = ClimaAnalysis.select(
                        era5_var_t;
                        by = ClimaAnalysis.MatchValue(),
                        pressure,
                    )
                    # Sometimes the float type of the dims don't match so we resample...
                    # sim_var_t_p = ClimaAnalysis.resampled_as(sim_var_t_p, era5_var_t_p)
                    ClimaAnalysis.Visualize.plot_bias_on_globe!(
                        fig[i, j],
                        sim_var_t_p,
                        era5_var_t_p,
                        # era5_var_t_p;
                        # cmap_extrema,
                    )
                end
            else
                ClimaAnalysis.Visualize.plot_bias_on_globe!(
                    fig[i, 1],
                    sim_var_t,
                    era5_var_t;
                    cmap_extrema,
                )
            end
        catch e
            @error "bias plot error: $(ClimaAnalysis.short_name(sim_var_t))"
        end
    end

    GeoMakie.save(joinpath(output_dir, "bias_sample_dates.png"), fig)
    return nothing
end

"""
    ClimaCalibrate.analyze_iteration(
        interface::CouplerModelInterface,
        ekp,
        g_ensemble,
        prior,
        output_dir,
        iteration,
    )

Analyze each iteration is completed by
- plotting the prior,
- plotting the contrained parameters and errors,
- plotting the G ensemble against the observation,
- reporting the residual diagnostics,
- plotting the bias,
- computing the ensemble spread.

Every plot is guarded so that one failing plot does not take the others with it.
ClimaCalibrate already catches a failure of this function, so the calibration continues
either way.
"""
function ClimaCalibrate.analyze_iteration(
    interface::CouplerModelInterface,
    ekp,
    g_ensemble,
    prior,
    output_dir,
    iteration,
)
    plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
    guard(f, what) =
        try
            f()
        catch e
            @error "$what failed" exception = (e, catch_backtrace())
        end

    guard(() -> plot_prior(output_dir, prior), "Prior plotting")
    guard(
        () -> plot_constrained_params_and_errors(output_dir, ekp, prior),
        "Parameter and error plotting",
    )
    guard(() -> plot_g_ensemble(plot_output_path, ekp, iteration), "G ensemble plotting")
    guard(() -> report_residual(plot_output_path, ekp, iteration), "Residual diagnostics")

    (; config) = interface
    job_id = get_job_id(config)
    member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1)
    simdir_path = joinpath(member_path, job_id, "output_active")
    guard("Bias plotting") do
        simdir = ClimaAnalysis.SimDir(simdir_path)
        plot_bias_weekly(ekp, simdir, iteration; output_dir = plot_output_path)
    end

    @info "Ensemble spread: $(scalar_spread(ekp))"
    return nothing
end

"""
    plot_prior(output_dir, prior)

Plot the marginal distribution of every prior dimension and save it to `output_dir`.

The prior does not change between iterations, so this overwrites the same file. It is
worth having next to the posterior plots: a posterior pressed against the edge of its
prior means the prior, not the data, is setting the answer.
"""
function plot_prior(output_dir, prior)
    fig = CairoMakie.Figure(size = (800, 600))
    EKP.Visualize.plot_parameter_distribution(fig[1, 1], prior)
    CairoMakie.save(joinpath(output_dir, "prior.png"), fig)
    return nothing
end

"""
    plot_g_ensemble(output_dir, ekp, iteration)

Plot the ensemble members of the forward map, their mean, and the observation for
`iteration`, and save it to `output_dir`.

This is the plot that says whether the ensemble is approaching the observation at all, and
where along the observation vector it is not.
"""
function plot_g_ensemble(output_dir, ekp, iteration)
    fig = CairoMakie.Figure(size = (1000, 400))
    ax = CairoMakie.Axis(
        fig[1, 1],
        title = "Forward map and observation, iteration $iteration",
        xlabel = "Index in the observation vector",
        ylabel = "Value",
    )
    g = ClimaCalibrate.Visualization.plot_g!(
        ax,
        ekp;
        iter = iteration,
        color = :black,
        alpha = 0.2,
    )
    g_mean =
        ClimaCalibrate.Visualization.plot_g_mean!(ax, ekp; iter = iteration, color = :black)
    obs = ClimaCalibrate.Visualization.plot_obs!(ax, ekp; iter = iteration, color = :blue)
    CairoMakie.Legend(fig[1, 2], [g, g_mean, obs], ["members", "mean", "observation"])
    CairoMakie.save(joinpath(output_dir, "g_ensemble.png"), fig)
    return nothing
end

"""
    report_residual(output_dir, ekp, iteration)

Report how much of the remaining residual is structured rather than noise-like, and save
a per-variable bar chart to `output_dir`.

`analyze_residual` projects `obs - mean(G)` onto the leading eigenvectors of the noise
covariance and normalises by the eigenvalues, so the projections are z-scores. A
structured energy near one is what the noise model predicts; much larger than one means
the residual has structure the noise model does not account for, and the per-variable
split says which observation or which part of the observation map to look at.
"""
function report_residual(output_dir, ekp, iteration; n_eigenvectors = 3)
    result = ClimaCalibrate.analyze_residual(ekp, iteration; n_eigenvectors)
    # `result.metadata` is one ClimaAnalysis Metadata per variable, in the same order as
    # the per-variable vectors below.
    names = String[ClimaAnalysis.short_name(m) for m in result.metadata]
    @info "Residual diagnostics, iteration $iteration" structured_energy =
        result.structured_energy
    for (name, energy, norm) in
        zip(names, result.structured_energy_by_variable, result.residual_norm_by_variable)
        @info "  $name" structured_energy = energy residual_norm = norm
    end

    fig = CairoMakie.Figure(size = (900, 400))
    ax1 = CairoMakie.Axis(
        fig[1, 1],
        title = "Structured energy by variable (1 = noise model)",
        xticks = (1:length(names), names),
        xticklabelrotation = pi / 4,
    )
    CairoMakie.barplot!(ax1, 1:length(names), collect(result.structured_energy_by_variable))
    CairoMakie.hlines!(ax1, [1.0]; color = :red, linestyle = :dash)
    ax2 = CairoMakie.Axis(
        fig[1, 2],
        title = "Residual norm by variable",
        xticks = (1:length(names), names),
        xticklabelrotation = pi / 4,
    )
    CairoMakie.barplot!(ax2, 1:length(names), collect(result.residual_norm_by_variable))
    CairoMakie.save(joinpath(output_dir, "residual_diagnostics.png"), fig)
    return nothing
end

"""
    plot_constrained_params_and_errors(output_dir, ekp, prior)

Plot the constrained parameters and errors from `ekp` and `prior` and save
them to `output_dir`.
"""
function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 1], ekp, error_metric = "loss")
    EKP.Visualize.plot_error_over_time(fig[1, dim_size + 2], ekp, error_metric = "loss")
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
end

"""
    scalar_spread(ekp)

Compute the mean over ensemble members of the squared Euclidean distance of the
forward model outputs from the ensemble mean.
"""
function scalar_spread(ekp)
    g_mean_final = EKP.get_g_mean_final(ekp)
    g_final = EKP.get_g_final(ekp)
    sq_dists = [sum((col .- g_mean_final) .^ 2) for col in eachcol(g_final)]
    return Statistics.mean(sq_dists)
end
