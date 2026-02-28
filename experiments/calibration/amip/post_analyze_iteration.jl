include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "subseasonal",
        "observation_map.jl",
    ),
)

# Extend bias_plot_extrema (defined in subseasonal/observation_map.jl) with weekly variables
bias_plot_extrema["ta_850hPa"] = (-2, 2)
bias_plot_extrema["ta_500hPa"] = (-2, 2)
bias_plot_extrema["ta_200hPa"] = (-2, 2)
bias_plot_extrema["hur_850hPa"] = (-2, 2)
bias_plot_extrema["hur_500hPa"] = (-2, 2)
bias_plot_extrema["hur_200hPa"] = (-2, 2)

# TODO: Unify this with `plot_bias` in the subseasonal observation_map.jl
"""
    plot_bias_weekly(ekp, simdir, iteration; output_dir)

Plot bias maps comparing simulation output to ERA5 observations for all variables in
`CALIBRATE_CONFIG.short_names`. ERA5 vars are reconstructed from the EKP observation
object and denormalized to physical units when normalization is enabled.
"""
function plot_bias_weekly(ekp, simdir, iteration; output_dir = simdir.simulation_path)
    sample_date_range = CALIBRATE_CONFIG.sample_date_ranges[iteration + 1]
    calib_start, _ = sample_date_range

    # Reconstruct ERA5 OutputVars from the EKP observation object
    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs = ClimaCalibrate.ObservationRecipe.get_observations_for_nth_iteration(
        obs_series,
        iteration + 1,
    )

    # ERA5 vars are kept as-is (normalized), matching what enters the loss function
    era5_vars =
        mapreduce(ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, minibatch_obs)

    sim_vars = map(CALIBRATE_CONFIG.short_names) do short_name
        var = preprocess_var(get_var(short_name, simdir), sample_date_range)
    end

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

    fig = GeoMakie.Figure(size = (1500, 500 * length(var_pairs)))
    for (i, (sim_var, era5_var)) in enumerate(var_pairs)
        sn = ClimaAnalysis.short_name(sim_var)
        sim_var_t = ClimaAnalysis.select(sim_var; by = MatchValue(), time = calib_start)
        era5_var_t = ClimaAnalysis.select(era5_var; by = MatchValue(), time = calib_start)
        cmap_extrema = get(bias_plot_extrema, sn, extrema(sim_var_t.data))
        if ClimaAnalysis.has_pressure(sim_var_t)
            for (j, pressure) in enumerate(ClimaAnalysis.pressures(sim_var_t))
                sim_var_t_p = ClimaAnalysis.select(sim_var_t; by = MatchValue(), pressure)
                era5_var_t_p = ClimaAnalysis.select(era5_var_t; by = MatchValue(), pressure)
                # Sometimes the float type of the dims don't match so we resample...
                # sim_var_t_p = ClimaAnalysis.resampled_as(sim_var_t_p, era5_var_t_p)
                ClimaAnalysis.Visualize.plot_bias_on_globe!(
                    fig[i, j],
                    sim_var_t_p,
                    era5_var_t_p;
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
    end

    GeoMakie.save(joinpath(output_dir, "bias_sample_dates.png"), fig)
    return nothing
end

function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
    plot_constrained_params_and_errors(output_dir, ekp, prior)

    job_id = get_job_id()
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
