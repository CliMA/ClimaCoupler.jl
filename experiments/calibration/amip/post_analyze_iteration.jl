bias_plot_extrema = Dict(
    "tas" => (-6, 6),
    "tas - ta" => (-6, 6),
    "hfls" => (-50, 50),
    "hfss" => (-25, 25),
    "rsus" => (-50, 50),
    "rlus" => (-50, 50),
    "mslp" => (-1000, 1000),
    # `pr` is compared in mm/day (the GPCP unit; sim is converted from
    # kg m^-2 s^-1 in preprocess_sim_vars), so these bounds are mm/day. The old
    # (-1e-4, 1e-4) was for the raw kg m^-2 s^-1 field and is off by ~86400.
    "pr" => (-3, 3),
    "lwp" => (-0.05, 0.05),
    "cl" => (-0.3, 0.3),
    "ta_850hPa" => (-2, 2),
    "ta_500hPa" => (-2, 2),
    "ta_200hPa" => (-2, 2),
    "hur_850hPa" => (-2, 2),
    "hur_500hPa" => (-2, 2),
    "hur_200hPa" => (-2, 2),
)

"""
    level_dim_of(var)

Return `(name, values)` for the vertical coordinate of `var` — `("pressure", ...)`,
`("altitude", ...)`, or `(nothing, [nothing])` when the variable has no vertical
dimension. Lets the bias plotting loop over levels without caring which kind of
vertical coordinate a variable uses (`ta`/`hur` are on pressure levels, `cl` is on
altitude levels, `lwp`/`pr` are 2-D).
"""
function level_dim_of(var)
    ClimaAnalysis.has_pressure(var) && return ("pressure", ClimaAnalysis.pressures(var))
    ClimaAnalysis.has_altitude(var) && return ("altitude", ClimaAnalysis.altitudes(var))
    return (nothing, [nothing])
end

"""
    select_level(var, level_name, value)

Select a single vertical level from `var`, or return `var` unchanged when
`level_name` is `nothing`.
"""
function select_level(var, level_name, value)
    isnothing(level_name) && return var
    return ClimaAnalysis.select(
        var;
        by = ClimaAnalysis.MatchValue(),
        NamedTuple{(Symbol(level_name),)}((value,))...,
    )
end

"""
    latitude_profile(var)

Reduce `var` to `(latitudes, values)` with `values` a plain vector aligned to
`latitudes`, squeezing any remaining singleton dimensions (e.g. a selected time or
level). Used for the zonal-mean bias panels, where the variable has no longitude.
"""
function latitude_profile(var)
    lats = ClimaAnalysis.latitudes(var)
    lat_idx = var.dim2index[ClimaAnalysis.latitude_name(var)]
    data = Array(var.data)
    # Move latitude first, then flatten the rest (all singleton after selection).
    perm = [lat_idx; setdiff(1:ndims(data), lat_idx)]
    permuted = reshape(permutedims(data, perm), size(data, lat_idx), :)
    return lats, permuted[:, 1]
end

"""
    plot_zonal_bias!(gridpos, sim_var, obs_var; title)

Plot a zonal-mean bias panel at `gridpos`: the simulation and observation as
latitude profiles (top), and their difference `sim - obs` (bottom) with a zero
line.

This is the no-longitude counterpart to `plot_bias_on_globe!`. The calibration
observations are zonal means (see `zonal_average` in preprocessing.jl), so they
have no longitude dimension and cannot be drawn on a globe at all — which is why
the on-globe path had been failing for every variable in every iteration. Latitude
is the only remaining spatial coordinate, so the bias *is* a latitude profile.

The reported mean bias and RMSE are area-weighted by `cos(latitude)`, since equal
latitude bands do not represent equal area; an unweighted mean would let the poles
dominate.
"""
function plot_zonal_bias!(gridpos, sim_var, obs_var; title)
    sim_lats, sim_vals = latitude_profile(sim_var)
    obs_lats, obs_vals = latitude_profile(obs_var)
    sim_lats == obs_lats || length(sim_lats) == length(obs_lats) || error(
        "simulation and observation latitude grids differ ($(length(sim_lats)) vs " *
        "$(length(obs_lats))); they must be on a common grid to take a bias",
    )

    bias = sim_vals .- obs_vals
    weights = cosd.(sim_lats)
    keep = findall(i -> isfinite(bias[i]), eachindex(bias))
    if isempty(keep)
        @warn "no overlapping finite points for zonal bias" title
        return nothing
    end
    wsum = sum(weights[keep])
    mean_bias = sum(bias[keep] .* weights[keep]) / wsum
    rmse = sqrt(sum(abs2.(bias[keep]) .* weights[keep]) / wsum)

    units = ClimaAnalysis.units(sim_var)
    gl = CairoMakie.GridLayout(gridpos)
    ax1 = CairoMakie.Axis(
        gl[1, 1],
        title = "$title\narea-weighted bias = $(round(mean_bias; sigdigits = 3)), " *
                "RMSE = $(round(rmse; sigdigits = 3)) [$units]",
        ylabel = "value [$units]",
    )
    sim_line = CairoMakie.lines!(ax1, sim_lats, sim_vals; color = :black)
    obs_line = CairoMakie.lines!(ax1, obs_lats, obs_vals; color = :blue)
    CairoMakie.axislegend(
        ax1,
        [sim_line, obs_line],
        ["simulation", "observation"];
        position = :lt,
    )
    CairoMakie.hidexdecorations!(ax1; grid = false)

    ax2 = CairoMakie.Axis(
        gl[2, 1],
        xlabel = "latitude [deg]",
        ylabel = "bias [$units]",
    )
    CairoMakie.hlines!(ax2, [0.0]; color = :black)
    CairoMakie.lines!(ax2, sim_lats, bias; color = :red)
    # Deliberately NOT clamped to `bias_plot_extrema`: those bounds are chosen to
    # keep an on-globe colour scale readable, and imposing them here clipped the
    # largest biases (for lwp the ±0.05 bound cut off both the tropical spike and
    # the high-latitude deficit) — exactly the points worth seeing. Let the axis
    # autoscale to the data instead.
    CairoMakie.rowsize!(gl, 1, CairoMakie.Relative(0.6))
    CairoMakie.linkxaxes!(ax1, ax2)
    return nothing
end

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

    # One row per variable, one column per vertical level (a single column for 2-D
    # variables like lwp/pr).
    n_cols = maximum(length(last(level_dim_of(sim_var))) for (sim_var, _) in var_pairs)
    fig = GeoMakie.Figure(size = (1000 * n_cols, 550 * length(var_pairs)))
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
        cmap_extrema = get(bias_plot_extrema, sn, nothing)
        # `lwp` (MAC) is an ocean-only retrieval with NaNs over land. The on-globe
        # plot's internal resampling is not NaN-aware, so mask to ocean to align
        # the NaN patterns. Only meaningful for the on-globe path.
        plot_mask = sn == "lwp" ? ClimaAnalysis.Visualize.oceanmask() : nothing

        # Loop the vertical coordinate generically: `pressure` for ta/hur,
        # `altitude` for cl, and a single no-op level for 2-D fields. Previously
        # only `pressure` was handled, so an altitude-resolved variable such as cl
        # fell through to the 2-D branch and failed.
        level_name, level_values = level_dim_of(sim_var_t)
        for (j, level) in enumerate(level_values)
            label = isnothing(level_name) ? sn :
                    "$sn @ $(round(level; sigdigits = 4)) $level_name"
            try
                sim_l = select_level(sim_var_t, level_name, level)
                obs_l = select_level(era5_var_t, level_name, level)
                # Dispatch on whether a longitude dimension survived preprocessing.
                # The calibration observations are zonal means, so in practice this
                # takes the latitude-profile path; the on-globe path is kept for
                # configurations that skip zonal averaging.
                if ClimaAnalysis.has_longitude(sim_l) &&
                   ClimaAnalysis.has_longitude(obs_l)
                    ClimaAnalysis.Visualize.plot_bias_on_globe!(
                        fig[i, j],
                        sim_l,
                        obs_l;
                        cmap_extrema = something(cmap_extrema, extrema(sim_l.data)),
                        mask = plot_mask,
                    )
                else
                    plot_zonal_bias!(fig[i, j], sim_l, obs_l; title = label)
                end
            catch e
                @error "bias plot error: $label" exception = (e, catch_backtrace())
            end
        end
    end

    figpath = joinpath(output_dir, "bias_sample_dates.png")
    GeoMakie.save(figpath, fig)
    @info "Saved bias plot" figpath
    return nothing
end

"""
    plot_g_vs_obs(ekp, iteration; g_ensemble = nothing, output_dir)

Plot the forward-model output against the observation POINT BY POINT in the
flattened observation index space, and save `g_vs_obs.png`.

Two rows:
 1. the G ensemble members, the ensemble mean forward map, and the observation,
    using the `ClimaCalibrate.Visualization` recipes;
 2. the residual `(mean(G) - y) / σ`, normalized by each point's observational
    noise standard deviation (the diagonal of the observation covariance), with
    ±1σ/±2σ guides and dashed lines at the per-variable block boundaries.

Why this exists: `plot_bias_weekly` draws bias on the globe, which requires a
longitude dimension — but the calibration observations are ZONAL MEANS, so that
plot fails for every variable (it has been erroring as "bias plot error: lwp" in
every iteration and saving blank panels). This diagnostic always works regardless
of the observation's dimensionality, since it operates on the flattened vector EKP
actually sees.

Row 1 uses the `ClimaCalibrate.Visualization.plot_g!` / `plot_g_mean!` /
`plot_obs!` recipes. Those live in `ClimaCalibrateMakieExt`, which the AMIP
Manifest had failed to register as an extension of the pinned `ne/async`
ClimaCalibrate even though the package declares it (and Makie 0.24 / CairoMakie
0.15 already satisfy its compat) — the Manifest entry was stale, so the recipe
methods silently never loaded. Registering the extension there makes them
available.

Note the recipes read the G ensemble as `EKP.get_g(ekp, iter)`, so `iter` must be
an iteration the ekp has actually completed. `analyze_iteration` runs *after*
`update_ensemble!`, so live `iter = iteration` is correct. A saved
`iteration_N/eki_file.jld2`, by contrast, is the state BEFORE iteration N runs and
holds only N-1 G ensembles, so re-plotting from one must clamp — see `g_iter`.

The normalized residual is the useful quantity, not the raw one: the observation
vector concatenates variables in different physical units (lwp ~0.09 kg m^-2 next
to pr ~-2.3 mm/day), so a raw residual is dominated by whichever variable has the
larger magnitude. Dividing by σ puts every point on a common "how many noise
standard deviations off am I" scale, which is exactly the quantity the EKP
objective is built from — points beyond ±2σ are real, learnable misfit, while
scatter within ±1σ is noise the calibration should not chase.
"""
function plot_g_vs_obs(ekp, iteration; g_ensemble = nothing, output_dir)
    # Which iteration's G are we plotting? Normally `analyze_iteration` hands us
    # this iteration's G directly. When called without one (e.g. re-plotting from
    # a saved eki_file), fall back to the newest G the ekp actually holds: a saved
    # `iteration_N/eki_file.jld2` is the state BEFORE iteration N runs, so it
    # contains only N-1 completed G ensembles and asking for N would throw.
    # G and the observation must come from the SAME iteration, otherwise a
    # multi-year minibatch would compare a model year against another year's obs.
    g_iter = iteration
    g_ens = g_ensemble
    if isnothing(g_ens)
        n_completed = EKP.get_N_iterations(ekp)
        n_completed >= 1 ||
            error("ekp holds no completed G ensemble; nothing to plot")
        g_iter = min(iteration, n_completed)
        g_ens = EKP.get_g(ekp, g_iter)
    end

    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs =
        ClimaCalibrate.get_observations_for_nth_iteration(obs_series, g_iter)
    y = mapreduce(EKP.get_obs, vcat, minibatch_obs)

    # Per-point noise sigma = sqrt of the observation covariance diagonal, read
    # entry by entry so we never materialize a large dense covariance.
    sigma = mapreduce(vcat, minibatch_obs) do obs
        cov = EKP.get_obs_noise_cov(obs)
        [sqrt(abs(cov[i, i])) for i in 1:size(cov, 1)]
    end

    # Variable names and block lengths, so the residual panel can be read per
    # variable. Lengths use the same NaN-dropping flatten as the observation.
    obs_vars =
        mapreduce(ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, minibatch_obs)
    var_names = ClimaAnalysis.short_name.(obs_vars)
    var_lengths = [length(ClimaAnalysis.flatten(v).data) for v in obs_vars]

    # NaN-aware ensemble mean, so one failed member does not blank the curve.
    g_mean = map(eachrow(g_ens)) do row
        finite = filter(isfinite, row)
        isempty(finite) ? NaN : Statistics.mean(finite)
    end

    fig = CairoMakie.Figure(size = (1400, 900))
    ax1 = CairoMakie.Axis(
        fig[1, 1],
        title = "G ensemble, mean forward map, and observations (iteration $g_iter)",
        xlabel = "Index",
        ylabel = "Value",
    )
    g_plot = ClimaCalibrate.Visualization.plot_g!(
        ax1,
        ekp;
        iter = g_iter,
        color = :black,
        alpha = 0.2,
    )
    g_mean_plot = ClimaCalibrate.Visualization.plot_g_mean!(
        ax1,
        ekp;
        iter = g_iter,
        color = :black,
    )
    obs_plot =
        ClimaCalibrate.Visualization.plot_obs!(ax1, ekp; iter = g_iter, color = :blue)
    CairoMakie.Legend(
        fig[1, 2],
        [g_plot, g_mean_plot, obs_plot],
        ["G", "G mean", "Observation"],
    )
    # Mark where one variable's block ends and the next begins. Row 1 shares a
    # single y-axis across variables in different units (lwp ~0.09 kg m^-2 vs
    # pr ~-2.3 mm/day), so the smaller-magnitude variable is squashed there —
    # read row 2 for anything cross-variable.
    let off = 0
        for len in var_lengths[1:(end - 1)]
            off += len
            CairoMakie.vlines!(ax1, [off + 0.5]; color = :black, linestyle = :dash)
        end
    end

    n = min(length(g_mean), length(y), length(sigma))
    resid = [
        (isfinite(g_mean[i]) && sigma[i] > 0) ? (g_mean[i] - y[i]) / sigma[i] : NaN
        for i in 1:n
    ]

    ax2 = CairoMakie.Axis(
        fig[2, 1],
        title = "Point-by-point residual (mean(G) - obs) / σ",
        xlabel = "Index",
        ylabel = "Residual [σ]",
    )
    CairoMakie.hlines!(ax2, [0.0]; color = :blue)
    CairoMakie.hlines!(ax2, [-1.0, 1.0]; color = :gray, linestyle = :dash)
    CairoMakie.hlines!(ax2, [-2.0, 2.0]; color = :red, linestyle = :dot)
    CairoMakie.scatter!(ax2, 1:n, resid; color = :black, markersize = 5)

    # Mark and label per-variable blocks, and report each block's RMS residual
    # in sigma units -- ~1 means "already fit to within noise" (nothing to
    # learn), >>1 means real remaining signal.
    offset = 0
    labels = String[]
    for (name, len) in zip(var_names, var_lengths)
        lo, hi = offset + 1, min(offset + len, n)
        offset += len
        offset < n && CairoMakie.vlines!(ax2, [offset + 0.5]; color = :black, linestyle = :dash)
        block = filter(isfinite, view(resid, lo:hi))
        rms = isempty(block) ? NaN : sqrt(Statistics.mean(abs2, block))
        push!(labels, "$name: RMS = $(round(rms; digits = 2))σ (n=$(hi - lo + 1))")
    end
    ax2.subtitle = join(labels, "   |   ")

    figpath = joinpath(output_dir, "g_vs_obs.png")
    CairoMakie.save(figpath, fig)
    @info "Saved point-by-point residual plot" figpath residual_rms_per_var = labels
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

    # Point-by-point residual in the flattened observation space. Kept separate
    # from (and ahead of) the on-globe bias plot below because that one requires a
    # longitude dimension and therefore fails on our zonal-mean observations —
    # this is the residual diagnostic we always get.
    try
        plot_g_vs_obs(ekp, iteration; g_ensemble, output_dir = plot_output_path)
    catch e
        @error "G-vs-obs residual plotting failed" exception = (e, catch_backtrace())
    end

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
