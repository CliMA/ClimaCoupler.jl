import ClimaAnalysis
import Dates
import ClimaCoupler
import CairoMakie
import GeoMakie
include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/data_sources.jl"))

"""
    get_cre_var(output_dir)
"""
function get_sim_cre_var(iteration_dir)
    # Need to recursively go through each iteration and get the diagnostics
    # corresponding to rsut and rsutcs
    member_paths = readdir(iteration_dir, join = true)
    filter!(x -> occursin("member_001", x), member_paths)
    member_paths = joinpath.(member_paths, "model_config", "output_active", "clima_atmos")

    # Code below is for taking average over all the iterations
    # simdirs = ClimaAnalysis.SimDir.(member_paths)
    # rsut_vars = [get(simdir; short_name = "rsut", reduction = "average", period = "1M") for simdir in simdirs]
    # rsutcs_vars = [get(simdir; short_name = "rsutcs", reduction = "average", period = "1M") for simdir in simdirs]
    # @assert length(rsut_vars) == length(rsutcs_vars)
    # cre_vars = [rsutcs_var - rsut_var for (rsutcs_var, rsut_var) in zip(rsutcs_vars, rsut_vars)]
    # mean_cre_data = sum(cre_var.data for cre_var in cre_vars) ./ length(cre_vars)
    # cre_var = ClimaAnalysis.remake(first(cre_vars), data = mean_cre_data)
    # cre_var = ClimaAnalysis.shift_to_start_of_previous_month(cre_var)
    # cre_var = ClimaAnalysis.set_units(cre_var, "W m^-2")
    # return cre_var

    simdir = ClimaAnalysis.SimDir(first(member_paths))
    sim_var_rsut = get(simdir, short_name = "rsut")
    sim_var_rsutcs = get(simdir, short_name = "rsutcs")
    sim_var_cre = sim_var_rsutcs - sim_var_rsut
    sim_var_cre = ClimaAnalysis.shift_to_start_of_previous_month(sim_var_cre)
    sim_var_cre = ClimaAnalysis.set_units(sim_var_cre, "W m^-2")
    return sim_var_cre
end

function get_obs_cre_var(start_date)
    rad_and_pr_obs_dict = get_obs_var_dict()
    rsut = rad_and_pr_obs_dict["rsut"](start_date)
    rsutcs = rad_and_pr_obs_dict["rsutcs"](start_date)
    cre = rsutcs - rsut
    cre = ClimaAnalysis.set_units(cre, "W m^-2")
    return cre
end

"""
    plot_cre_leaderboard(var, spinup)

Plot observational data, simulation data, and bias.
"""
function plot_cre!(fig, iteration_dir, spinup)
    # Get sim and obs
    sim_var = get_sim_cre_var(iteration_dir)
    start_date = Dates.DateTime(sim_var.attributes["start_date"])
    month_in_secs = 32 * 86_400.0

    sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup * month_in_secs)
    obs_var = get_obs_cre_var(start_date)
    obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)
    sim_var = ClimaAnalysis.average_time(sim_var)
    obs_var = ClimaAnalysis.average_time(obs_var)

    # TODO: Add color range here...
    ClimaAnalysis.Visualize.plot_bias_on_globe!(fig, sim_var, obs_var)
    supertitle = CairoMakie.Label(fig[0, :], last(splitpath(iteration_dir)), fontsize = 30, tellheight = true)
end

"""
    plot_cre_ann_seasons!(fig, iteration_dir, spinup; season_names = false)

Plot a row of annual and seasonal biases.

If `season_names` is `true`, then the row above will include the names of the
seasons.
"""
function plot_cre_ann_seasons!(fig, iteration_dir, spinup; plot_season_names = false)
    # Get sim and obs
    sim_var = get_sim_cre_var(iteration_dir)
    start_date = Dates.DateTime(sim_var.attributes["start_date"])
    month_in_secs = 32 * 86_400.0

    sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup * month_in_secs)
    obs_var = get_obs_cre_var(start_date)
    obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

    avg_time(season) = isempty(season) ? season : ClimaAnalysis.average_time(season)
    obs_seasons = avg_time.(ClimaAnalysis.split_by_season(obs_var))
    sim_seasons = avg_time.(ClimaAnalysis.split_by_season(sim_var))

    sim_var = ClimaAnalysis.average_time(sim_var)
    obs_var = ClimaAnalysis.average_time(obs_var)

    obs_ann_seasons = (obs_var, obs_seasons...)
    sim_ann_seasons = (sim_var, sim_seasons...)
    if plot_season_names
        season_names = ("ANN", (isempty(var) ? "" : var.attributes["season"] for var in sim_seasons)...)
        season_names = filter(x -> x != "", season_names)
    end
    i = 1
    for (sim, obs) in zip(sim_ann_seasons, obs_ann_seasons)
        isempty(sim) && continue
        layout = fig[1, i] = CairoMakie.GridLayout()
        ClimaAnalysis.Visualize.plot_bias_on_globe!(layout, sim, obs)
        if plot_season_names
            CairoMakie.Label(layout[0, 1], season_names[i], fontsize = 30, tellwidth = false, justification = :center)
        end
        i += 1
    end
    CairoMakie.Label(fig[1, 0], last(splitpath(iteration_dir)), fontsize = 30, tellheight = false)
    return nothing
end

"""
    find_ANN_and_season_names(iteration_dir, spinup)

Find the seasons that is present in the simulation data.
"""
function find_ANN_and_season_names(iteration_dir, spinup)
    sim_var = get_sim_cre_var(iteration_dir)
    month_in_secs = 32 * 86_400.0
    sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup * month_in_secs)

    sim_seasons = ClimaAnalysis.split_by_season(sim_var)
    season_names = ("ANN", (isempty(var) ? "" : var.attributes["season"] for var in sim_seasons)...)
    return filter(x -> x != "", season_names)
end

"""
    plot_cre_leaderboard_from_iters(output_dir, spinup, iters)

Plot leaderboard for each iteration in a square and a seasonal leaderboard.
"""
function plot_cre_leaderboard_from_iters(output_dir, spinup, iters)
    @info "Plotting CRE leaderboard from $iters iterations in $output_dir"
    iter_paths = readdir(output_dir, join = true)
    filter!(x -> occursin("iteration", x), iter_paths)
    iter_paths = iter_paths[1:iters]
    sort!(iter_paths)

    # Plot time average for each iteration
    n_iters = length(iter_paths)
    COLS = sqrt(n_iters) |> ceil
    fig = CairoMakie.Figure(size = (COLS * 600, COLS * 500))

    for (i, iter_path) in enumerate(iter_paths)
        r = Int((i - 1) ÷ COLS + 1)
        c = Int((i - 1) % COLS + 1)
        layout = fig[r, c] = CairoMakie.GridLayout()
        plot_cre!(layout, iter_path, spinup)
    end
    CairoMakie.save(joinpath(output_dir, "leaderboard.png"), fig)

    season_names = find_ANN_and_season_names(first(iter_paths), spinup)
    fig = CairoMakie.Figure(size = (800 * length(season_names), 500 * n_iters))
    for (i, iter_path) in enumerate(iter_paths)
        layout = fig[i, 1] = CairoMakie.GridLayout()
        if i == 1
            plot_cre_ann_seasons!(layout, iter_path, spinup, plot_season_names = true)
        else
            plot_cre_ann_seasons!(layout, iter_path, spinup)
        end
    end
    CairoMakie.save(joinpath(output_dir, "leaderboard_ann_seasons.png"), fig)
end
