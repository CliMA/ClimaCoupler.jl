include(joinpath(@__DIR__, "steering_indicators.jl"))
import Random
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
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
    "pr" => (-6, 6),
    "lwp" => (-0.1, 0.1),
    "cl" => (-0.3, 0.3),
    "clt" => (-0.5, 0.5),
    "swcre" => (-120, 120),
    "lwcre" => (-120, 120),
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
    snap_lonlat_coords(obs_var, sim_var)

Return `obs_var` with its longitude and latitude coordinate arrays replaced by
`sim_var`'s when the grids match in size and agree up to rounding. The two
grids are the same block means, but the observation stores them as Float32 and
the sim as Float64, so interpolation at the corner nodes can fall outside the
domain and error. Grids that differ by more than rounding are returned
unchanged.
"""
function snap_lonlat_coords(obs_var, sim_var)
    (ClimaAnalysis.has_longitude(obs_var) && ClimaAnalysis.has_longitude(sim_var)) ||
        return obs_var
    obs_lon = ClimaAnalysis.longitude_name(obs_var)
    obs_lat = ClimaAnalysis.latitude_name(obs_var)
    sim_lons = collect(Float64, ClimaAnalysis.longitudes(sim_var))
    sim_lats = collect(Float64, ClimaAnalysis.latitudes(sim_var))
    obs_lons = collect(Float64, obs_var.dims[obs_lon])
    obs_lats = collect(Float64, obs_var.dims[obs_lat])
    (length(obs_lons) == length(sim_lons) && length(obs_lats) == length(sim_lats)) ||
        return obs_var
    (isapprox(obs_lons, sim_lons; rtol = 1e-4) &&
     isapprox(obs_lats, sim_lats; rtol = 1e-4)) || return obs_var
    new_dims = copy(obs_var.dims)
    new_dims[obs_lon] = sim_lons
    new_dims[obs_lat] = sim_lats
    return ClimaAnalysis.remake(obs_var; dims = new_dims)
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
        # `lwp` (MAC) is an ocean-only retrieval with NaNs over land, so mask
        # out the continents. `landmask()` covers land; `oceanmask()` covers
        # the ocean and painted a uniform polygon over all the data. Only
        # meaningful for the on-globe path.
        plot_mask = sn == "lwp" ? ClimaAnalysis.Visualize.landmask() : nothing

        # Loop the vertical coordinate generically: `pressure` for ta/hur,
        # `altitude` for cl, and a single no-op level for 2-D fields. Previously
        # only `pressure` was handled, so an altitude-resolved variable such as cl
        # fell through to the 2-D branch and failed.
        level_name, level_values = level_dim_of(sim_var_t)
        for (j, level) in enumerate(level_values)
            label = isnothing(level_name) ? "$sn (iteration $iteration)" :
                    "$sn @ $(round(level; sigdigits = 4)) $level_name (iteration $iteration)"
            try
                sim_l = select_level(sim_var_t, level_name, level)
                obs_l = select_level(era5_var_t, level_name, level)
                # Dispatch on whether a longitude dimension survived preprocessing.
                # The calibration observations are zonal means, so in practice this
                # takes the latitude-profile path; the on-globe path is kept for
                # configurations that skip zonal averaging.
                if ClimaAnalysis.has_longitude(sim_l) &&
                   ClimaAnalysis.has_longitude(obs_l)
                    # Coarse-grid observations store Float32 coordinates while
                    # the coarsened sim carries Float64 block means. The corner
                    # nodes then differ by rounding, plot_bias_on_globe!'s
                    # internal resampling throws a BoundsError, and every
                    # panel is left empty. Snap the obs coordinates to the
                    # sim's before plotting.
                    obs_l = snap_lonlat_coords(obs_l, sim_l)
                    # plot_bias_on_globe!'s default more_kwargs sets
                    # :mask => Dict(), which drops the white fill and paints
                    # the mask polygon in Makie's default blue. Pass the
                    # full dict with an explicit white mask.
                    ClimaAnalysis.Visualize.plot_bias_on_globe!(
                        fig[i, j],
                        sim_l,
                        obs_l;
                        # The fallback must skip NaNs: coverage-masked fields
                        # (for example clt) contain them, and extrema over NaNs
                        # returns (NaN, NaN), which fails inside the plot call.
                        cmap_extrema = something(
                            cmap_extrema,
                            extrema(filter(isfinite, vec(sim_l.data))),
                        ),
                        mask = plot_mask,
                        more_kwargs = Dict(
                            :plot => Dict(),
                            :cb => Dict(),
                            :axis => Dict(),
                            :coast => Dict(:color => :black),
                            :mask => Dict(:color => :white),
                        ),
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

Plot the forward-model output against the observation point by point and save
`g_vs_obs.png`.

Each observed variable gets its own column, so variables with different scales
stay readable. Each column has two boxes:
 1. The G ensemble members, the ensemble mean, and the observation, on the
    variable's own value axis.
 2. The residual `(mean(G) - y) / σ` with guides at ±1σ and ±2σ, where σ is the
    square root of the observation covariance diagonal. The box title reports
    the block RMS in σ units. An RMS near 1 means the variable is fit to its
    noise level. An RMS well above 2 means learnable signal remains.

This plot works for observations of any dimensionality because it operates on
the flattened vector that EKP sees. `plot_bias_weekly` needs spatial
coordinates and is complementary.

`analyze_iteration` passes this iteration's G ensemble directly. When
`g_ensemble` is not given, the newest completed G in the ekp is used. A saved
`iteration_N/eki_file.jld2` holds the state from before iteration N runs, so it
contains only N-1 completed G ensembles. G and the observation always come from
the same iteration, so a minibatch never compares one sample's model output
against another sample's observation.
"""
function plot_g_vs_obs(ekp, iteration; g_ensemble = nothing, output_dir)
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

    # Read the covariance diagonal entry by entry to avoid materializing a
    # large dense matrix.
    sigma = mapreduce(vcat, minibatch_obs) do obs
        cov = EKP.get_obs_noise_cov(obs)
        [sqrt(abs(cov[i, i])) for i in 1:size(cov, 1)]
    end

    obs_vars =
        mapreduce(ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, minibatch_obs)
    var_names = ClimaAnalysis.short_name.(obs_vars)
    var_lengths = [length(ClimaAnalysis.flatten(v).data) for v in obs_vars]

    # Ensemble mean per point. NaN entries from failed members are skipped.
    g_mean = map(eachrow(g_ens)) do row
        finite = filter(isfinite, row)
        isempty(finite) ? NaN : Statistics.mean(finite)
    end

    n = min(length(g_mean), length(y), length(sigma))
    n_vars = length(var_names)
    fig = CairoMakie.Figure(size = (700 * n_vars, 900))

    labels = String[]
    legend_handles = nothing
    offset = 0
    for (col, (name, len)) in enumerate(zip(var_names, var_lengths))
        lo, hi = offset + 1, min(offset + len, n)
        offset += len
        lo <= hi || continue
        idx = lo:hi
        xs = 1:length(idx)

        # Box 1: values on this variable's own axis.
        ax1 = CairoMakie.Axis(
            fig[1, col],
            title = "$name (iteration $g_iter)",
            xlabel = "index in block",
            ylabel = "value",
        )
        h_g = nothing
        for member in eachcol(view(g_ens, idx, :))
            h_g = CairoMakie.lines!(ax1, xs, collect(member); color = (:black, 0.2))
        end
        h_mean = CairoMakie.lines!(ax1, xs, g_mean[idx]; color = :black)
        h_obs = CairoMakie.lines!(ax1, xs, y[idx]; color = :blue)
        isnothing(legend_handles) && (legend_handles = [h_g, h_mean, h_obs])

        # Box 2: residual in units of the observation noise.
        resid = [
            (isfinite(g_mean[i]) && sigma[i] > 0) ? (g_mean[i] - y[i]) / sigma[i] :
            NaN for i in idx
        ]
        finite_resid = filter(isfinite, resid)
        rms = isempty(finite_resid) ? NaN :
              sqrt(Statistics.mean(abs2, finite_resid))
        push!(labels, "$name: RMS = $(round(rms; digits = 2))σ (n=$(length(idx)))")
        ax2 = CairoMakie.Axis(
            fig[2, col],
            title = "residual (mean(G) - obs) / σ   RMS = $(round(rms; digits = 2))σ",
            xlabel = "index in block",
            ylabel = "residual [σ]",
        )
        CairoMakie.hlines!(ax2, [0.0]; color = :blue)
        CairoMakie.hlines!(ax2, [-1.0, 1.0]; color = :gray, linestyle = :dash)
        CairoMakie.hlines!(ax2, [-2.0, 2.0]; color = :red, linestyle = :dot)
        CairoMakie.scatter!(ax2, xs, resid; color = :black, markersize = 5)
    end

    isnothing(legend_handles) ||
        CairoMakie.Legend(fig[1, n_vars + 1], legend_handles,
                          ["G", "G mean", "Observation"])

    figpath = joinpath(output_dir, "g_vs_obs.png")
    CairoMakie.save(figpath, fig)
    @info "Saved point-by-point residual plot" figpath residual_rms_per_var = labels
    return nothing
end

"""
    plot_priors(priors, output_dir; n_samples = 100_000, rng_seed = 42)

Plot the marginal prior distribution of each parameter in constrained
(physical) space via EKP's Makie extension
(`EKP.Visualize.plot_parameter_distribution`) and save `priors.png` in
`output_dir`.
"""
function plot_priors(priors, output_dir; n_samples = 100_000, rng_seed = 42)
    n_params = length(PD.get_name(priors))
    n_cols = ceil(Int, sqrt(n_params))
    n_rows = cld(n_params, n_cols)
    fig = CairoMakie.Figure(size = (420 * n_cols, 320 * n_rows))
    EKP.Visualize.plot_parameter_distribution(
        fig,
        priors;
        constrained = true,
        n_sample = n_samples,
        rng = Random.MersenneTwister(rng_seed),
    )
    figpath = joinpath(output_dir, "priors.png")
    CairoMakie.save(figpath, fig)
    @info "Saved prior distributions plot" figpath
    return nothing
end

"""
    iteration_label(frame_path)

Column header for a plot living in an `iteration_NNN` directory:
"iteration N".
"""
function iteration_label(frame_path)
    itdir = basename(dirname(frame_path))
    n = tryparse(Int, last(split(itdir, '_')))
    return isnothing(n) ? itdir : "iteration $n"
end

"""
    animate_iteration_plots(output_dir; names, delay_cs = 100, strip = true)

Collect one plot type from every `iteration_NNN` directory and write an
animated GIF `<stem>_evolution.gif` in `output_dir`. With `strip = true`, also
write `<stem>_first_last.png`, the first and last frames side by side. Uses
ImageMagick, which is on the default PATH on Derecho. Call after the last
iteration, or standalone on any completed run directory:

    animate_iteration_plots("/path/to/output_dir")

`delay_cs` is the frame delay in centiseconds.
"""
function animate_iteration_plots(
    output_dir;
    names = ["bias_sample_dates.png", "g_vs_obs.png"],
    delay_cs = 100,
    strip = true,
)
    magick = Sys.which("magick")
    if isnothing(magick)
        @warn "ImageMagick not found on PATH; skipping iteration-plot animation"
        return nothing
    end
    for name in names
        frames = sort(
            filter(
                isfile,
                [
                    joinpath(d, name) for d in readdir(output_dir; join = true) if
                    startswith(basename(d), "iteration_")
                ],
            ),
        )
        if length(frames) < 2
            @info "Fewer than two frames for $name; skipping animation"
            continue
        end
        stem = splitext(name)[1]
        gif = joinpath(output_dir, "$(stem)_evolution.gif")
        try
            run(`$magick -delay $delay_cs $frames -loop 0 $gif`)
            @info "Wrote iteration animation" gif n_frames = length(frames)
            if strip
                strippath = joinpath(output_dir, "$(stem)_first_last.png")
                labeled(frame) = [
                    "(",
                    frame,
                    "-resize",
                    "x1600>",
                    "-background",
                    "white",
                    "-gravity",
                    "north",
                    "-splice",
                    "0x90",
                    "-pointsize",
                    "64",
                    "-fill",
                    "black",
                    "-annotate",
                    "+0+15",
                    iteration_label(frame),
                    ")",
                ]
                run(Cmd([
                    magick,
                    labeled(first(frames))...,
                    labeled(last(frames))...,
                    "+append",
                    strippath,
                ]))
                @info "Wrote first/last comparison" strippath
            end
        catch e
            @error "Failed to animate $name" exception = (e, catch_backtrace())
        end
    end
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

    # Log the plain-language steering block. Advisory and never fatal.
    try
        @info "\n" * steering_summary(ekp, g_ensemble, prior, iteration)
    catch e
        @error "Steering indicators failed" exception = (e, catch_backtrace())
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

    # Refresh the cross-iteration animations. The GIFs stay current even when a
    # run stops before its final iteration.
    try
        animate_iteration_plots(output_dir)
    catch e
        @error "Iteration-plot animation failed" exception = (e, catch_backtrace())
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
