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
        # `lwp` (MAC) is an ocean-only retrieval with NaNs over land. The bias
        # plot's internal resampling is not NaN-aware, which previously errored
        # ("bias plot error: lwp"). Masking the ocean aligns the NaN pattern of
        # both fields so the bias can be computed and plotted.
        plot_mask = sn == "lwp" ? ClimaAnalysis.Visualize.oceanmask() : nothing
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
                    mask = plot_mask,
                )
            end
        catch e
            @error "bias plot error: $(ClimaAnalysis.short_name(sim_var_t))" exception =
                (e, catch_backtrace())
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
- plotting the contrained parameters and errors,
- plotting the bias,
- computing the ensemble spread.
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
    plot_constrained_params_and_errors(output_dir, ekp, prior)

    (; config) = interface
    job_id = get_job_id(config)
    member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1)
    simdir_path = joinpath(member_path, job_id, "output_active")
    try
        simdir = ClimaAnalysis.SimDir(simdir_path)
        plot_bias_weekly(ekp, simdir, iteration; output_dir = plot_output_path)
    catch e
        @error "Bias plotting failed" exception = (e, catch_backtrace())
    end

    @info "Ensemble spread: $(scalar_spread(ekp))"
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
