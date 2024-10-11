import ClimaAnalysis
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import GeoMakie
import CairoMakie
import Dates

include("data_sources.jl")

"""
    compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)

Plot the biases and a leaderboard of various variables.

The paramter `leaderboard_base_path` is the path to save the leaderboards and bias plots and
`diagnostics_folder_path` is the path to the simulation data.
"""
function compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)
    @info "Error against observations"

    # Set up dict for storing simulation and observational data after processing
    sim_obs_comparsion_dict = Dict()
    seasons = ["ANN", "MAM", "JJA", "SON", "DJF"]

    # Print dates for debugging
    _, var_func = first(sim_var_dict)
    var = var_func()
    output_dates = Dates.DateTime(var.attributes["start_date"]) .+ Dates.Second.(ClimaAnalysis.times(var))
    @info "Working with dates:"
    @info output_dates

    for short_name in keys(sim_var_dict)
        # Simulation data
        sim_var = sim_var_dict[short_name]()

        # Observational data
        obs_var = obs_var_dict[short_name](sim_var.attributes["start_date"])

        # Remove first spin_up_months from simulation
        spin_up_months = 6
        spinup_cutoff = spin_up_months * 31 * 86400.0
        ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff &&
            (sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff))

        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)
        obs_var_seasons = ClimaAnalysis.split_by_season(obs_var)
        sim_var_seasons = ClimaAnalysis.split_by_season(sim_var)

        # Add annual to start of seasons
        obs_var_seasons = (obs_var, obs_var_seasons...)
        sim_var_seasons = (sim_var, sim_var_seasons...)

        # Take time average
        obs_var_seasons = (
            !ClimaAnalysis.isempty(obs_var) ? ClimaAnalysis.average_time(obs_var) : obs_var for
            obs_var in obs_var_seasons
        )
        sim_var_seasons = (
            !ClimaAnalysis.isempty(sim_var) ? ClimaAnalysis.average_time(sim_var) : sim_var for
            sim_var in sim_var_seasons
        )

        # Save observation and simulation data
        sim_obs_comparsion_dict[short_name] = Dict(
            season => (sim_var_s, obs_var_s) for
            (season, sim_var_s, obs_var_s) in zip(seasons, sim_var_seasons, obs_var_seasons)
        )
    end

    # Filter seasons to remove seasons with no dates
    _, var = first(sim_obs_comparsion_dict)
    filter!(season -> !ClimaAnalysis.isempty(var[season][1]), seasons)

    # Plot annual plots
    for compare_vars_biases in compare_vars_biases_groups
        fig_bias = CairoMakie.Figure(; size = (600, 300 * length(compare_vars_biases)))
        for (row, short_name) in enumerate(compare_vars_biases)
            sim_var = sim_obs_comparsion_dict[short_name]["ANN"][1]
            if !ClimaAnalysis.isempty(sim_var)
                # Add "mean " for plotting the title
                sim_var.attributes["short_name"] = "mean $(ClimaAnalysis.short_name(sim_var))"
                cmap_extrema = get(compare_vars_biases_plot_extrema, short_name) do
                    extrema(ClimaAnalysis.bias(sim_obs_comparsion_dict[short_name]["ANN"]...).data)
                end
                ClimaAnalysis.Visualize.plot_bias_on_globe!(
                    fig_bias,
                    sim_obs_comparsion_dict[short_name]["ANN"]...,
                    cmap_extrema = cmap_extrema,
                    p_loc = (row, 1),
                )
            end
        end
        CairoMakie.save(joinpath(leaderboard_base_path, "bias_$(first(compare_vars_biases))_ANN.png"), fig_bias)
    end

    # Plot plots with all the seasons (not annual)
    seasons_no_ann = filter(season -> season != "ANN", seasons)
    for compare_vars_biases in compare_vars_biases_groups
        fig_all_seasons = CairoMakie.Figure(; size = (600 * length(seasons_no_ann), 300 * length(compare_vars_biases)))
        for (col, season) in enumerate(seasons_no_ann)
            # Plot at 2 * col - 1 because a bias plot takes up 1 x 2 space
            CairoMakie.Label(fig_all_seasons[0, 2 * col - 1], season, tellwidth = false, fontsize = 30)
            for (row, short_name) in enumerate(compare_vars_biases)
                sim_var = sim_obs_comparsion_dict[short_name][season][1]
                if !ClimaAnalysis.isempty(sim_var)
                    cmap_extrema = get(compare_vars_biases_plot_extrema, short_name) do
                        extrema(ClimaAnalysis.bias(sim_obs_comparsion_dict[short_name][season]...).data)
                    end
                    # Add "mean " for plotting the title
                    sim_var.attributes["short_name"] = "mean $(ClimaAnalysis.short_name(sim_var))"
                    ClimaAnalysis.Visualize.plot_bias_on_globe!(
                        fig_all_seasons,
                        sim_obs_comparsion_dict[short_name][season]...,
                        cmap_extrema = cmap_extrema,
                        p_loc = (row, 2 * col - 1),
                    )
                end
            end
        end
        CairoMakie.save(
            joinpath(leaderboard_base_path, "bias_$(first(compare_vars_biases))_all_seasons.png"),
            fig_all_seasons,
        )
    end

    # Add RMSE for the CliMA model and for each season
    for short_name in rmse_var_names
        for season in seasons
            !ClimaAnalysis.isempty(sim_obs_comparsion_dict[short_name][season][1]) && (
                rmse_var_dict[short_name]["CliMA", season] =
                    ClimaAnalysis.global_rmse(sim_obs_comparsion_dict[short_name][season]...)
            )
        end
    end

    # Plot box plots
    fig_leaderboard = CairoMakie.Figure(; size = (800, 300 * length(rmse_var_dict) + 400), fontsize = 20)
    for (loc, short_name) in enumerate(rmse_var_names)
        ClimaAnalysis.Visualize.plot_boxplot!(
            fig_leaderboard,
            rmse_var_dict[short_name],
            ploc = (loc, 1),
            best_and_worst_category_name = "ANN",
        )
    end

    # Plot leaderboard
    ClimaAnalysis.Visualize.plot_leaderboard!(
        fig_leaderboard,
        (rmse_var_dict[short_name] for short_name in rmse_var_names)...,
        best_category_name = "ANN",
        ploc = (length(rmse_var_dict) + 1, 1),
    )
    CairoMakie.save(joinpath(leaderboard_base_path, "bias_leaderboard.png"), fig_leaderboard)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        error("Usage: julia leaderboard.jl <Filepath to save leaderboard and bias plots> <Filepath to simulation data>")
    end
    leaderboard_base_path = ARGS[begin]
    diagnostics_folder_path = ARGS[begin + 1]
    compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)
end
