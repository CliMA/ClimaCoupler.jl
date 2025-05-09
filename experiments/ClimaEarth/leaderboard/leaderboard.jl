import ClimaAnalysis
import GeoMakie
import CairoMakie
import Dates

include("data_sources.jl")

"""
    compute_leaderboard(leaderboard_base_path, diagnostics_folder_path, spinup)

Plot the biases and a leaderboard of various variables defined over longitude, latitude, and
time.

The argument `leaderboard_base_path` is the path to save the leaderboards and bias plots,
`diagnostics_folder_path` is the path to the simulation data, and `spinup` is the number
of months to remove from the beginning of the simulation.

Loading and preprocessing simulation data is done by `get_sim_var_dict`. Loading and
preprocessing observational data is done by `get_obs_var_dict`. The ranges of the bias plots
are determined by `get_compare_vars_biases_plot_extrema`. The groups of variables plotted on
the bias plots are determined by `get_compare_vars_biases_groups()`. Loading the RMSEs from
other models is done by `get_rmse_var_dict`. See the functions defined in data_sources.jl.
"""
function compute_leaderboard(leaderboard_base_path, diagnostics_folder_path, spinup)
    @info "Error against observations"

    # Get everything we need from data_sources.jl
    sim_var_dict = get_sim_var_dict(diagnostics_folder_path)
    obs_var_dict = get_obs_var_dict()
    compare_vars_biases_plot_extrema = get_compare_vars_biases_plot_extrema()
    rmse_var_dict = get_rmse_var_dict()
    compare_vars_biases_groups = get_compare_vars_biases_groups()

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
        spinup_months = spinup
        spinup_cutoff = spinup_months * 31 * 86400.0
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
        # Plot all the variables that we can
        compare_vars_biases = filter(x -> x in keys(sim_obs_comparsion_dict), compare_vars_biases)
        length(compare_vars_biases) > 0 || continue

        fig_bias = CairoMakie.Figure(; size = (600, 75.0 + 300 * length(compare_vars_biases)))
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
        # Plot all the variables that we can
        compare_vars_biases = filter(x -> x in keys(sim_obs_comparsion_dict), compare_vars_biases)
        length(compare_vars_biases) > 0 || continue

        fig_all_seasons =
            CairoMakie.Figure(; size = (600 * length(seasons_no_ann), 75.0 + 300 * length(compare_vars_biases)))
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
    is_in_sim_vars(k) = k in keys(sim_var_dict)
    rmse_var_dict = Dict(k => v for (k, v) in rmse_var_dict if is_in_sim_vars(k))
    rmse_var_names = keys(rmse_var_dict)
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

"""
    compute_pfull_leaderboard(leaderboard_base_path, diagnostics_folder_path, spinup)

Plot the biases and a leaderboard of various variables defined over longitude, latitude,
time, and pressure levels.

The argument `leaderboard_base_path` is the path to save the leaderboards and bias plots,
`diagnostics_folder_path` is the path to the simulation data, and `spinup` is the number
of months to remove from the beginning of the simulation.

Loading and preprocessing simulation data is done by `get_sim_var_in_pfull_dict`. Loading
and preprocessing observational data is done by `get_obs_var_in_pfull_dict`. The ranges of
the bias plots is defined by `get_compare_vars_biases_plot_extrema_pfull`. See the functions
defined in data_sources.jl for more information.
"""
function compute_pfull_leaderboard(leaderboard_base_path, diagnostics_folder_path, spinup)
    @info "Error against observations for variables in pressure coordinates"

    sim_var_pfull_dict = get_sim_var_in_pfull_dict(diagnostics_folder_path)
    obs_var_pfull_dict = get_obs_var_in_pfull_dict()

    # Print dates for debugging
    _, get_var = first(sim_var_pfull_dict)
    var = get_var()
    output_dates = Dates.DateTime(var.attributes["start_date"]) .+ Dates.Second.(ClimaAnalysis.times(var))
    @info "Working with dates:"
    @info output_dates

    # Set up dict for storing simulation and observational data after processing
    sim_obs_comparsion_dict = Dict()

    for short_name in keys(sim_var_pfull_dict)
        sim_var = sim_var_pfull_dict[short_name]()
        obs_var = obs_var_pfull_dict[short_name](sim_var.attributes["start_date"])

        # Check units for pressure are in hPa
        ClimaAnalysis.dim_units(sim_var, ClimaAnalysis.pressure_name(sim_var)) != "hPa" &&
            error("Units of pressure should be hPa for $short_name simulation data")
        ClimaAnalysis.dim_units(obs_var, ClimaAnalysis.pressure_name(obs_var)) != "hPa" &&
            error("Units of pressure should be hPa for $short_name simulation data")

        # Remove first spin_up_months from simulation
        spinup_cutoff = spinup * 31 * 86400.0
        ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff &&
            (sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff))

        # Restrain the pressure levels so we can resample
        min_pfull = ClimaAnalysis.pressures(obs_var)[1]
        sim_pressures = ClimaAnalysis.pressures(sim_var)
        greater_than_min_pfull_idx = findfirst(x -> x > min_pfull, sim_pressures)
        sim_var = ClimaAnalysis.window(sim_var, "pfull", left = sim_pressures[greater_than_min_pfull_idx])

        # Resample over lat, lon, time, and pressures
        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

        # Take time average
        sim_var = ClimaAnalysis.average_time(sim_var)
        obs_var = ClimaAnalysis.average_time(obs_var)

        # Save observation and simulation data
        sim_obs_comparsion_dict[short_name] = (; sim = sim_var, obs = obs_var)
    end

    # Slicing only uses the nearest value, so we also need to keep track of the
    # actual pressure levels that we get when slicing
    target_p_lvls = [850.0, 500.0, 250.0]
    real_p_lvls = []

    # Get units for pressure for plotting
    p_units = ClimaAnalysis.dim_units(first(sim_obs_comparsion_dict)[2].sim, "pfull")

    # Initialize ranges for colorbar and figure whose columns are pressure levels and rows are variables
    compare_vars_biases_plot_extrema_pfull = get_compare_vars_biases_plot_extrema_pfull()
    fig_bias_pfull_vars = CairoMakie.Figure(size = (650 * length(target_p_lvls), 450 * length(sim_obs_comparsion_dict)))

    # Plot bias at the pressure levels in p_lvls
    for (row_idx, short_name) in enumerate(keys(sim_obs_comparsion_dict))
        # Plot label for variable name
        CairoMakie.Label(fig_bias_pfull_vars[row_idx, 0], short_name, tellheight = false, fontsize = 30)
        for (col_idx, p_lvl) in enumerate(target_p_lvls)
            layout = fig_bias_pfull_vars[row_idx, col_idx] = CairoMakie.GridLayout()
            sim_var = sim_obs_comparsion_dict[short_name].sim
            obs_var = sim_obs_comparsion_dict[short_name].obs

            # Slice at pressure level using nearest value rather interpolating
            sim_var = ClimaAnalysis.slice(sim_var, pfull = p_lvl)
            obs_var = ClimaAnalysis.slice(obs_var, pressure_level = p_lvl)

            # Get the actual pressure level that we slice at
            push!(real_p_lvls, parse(Float64, sim_var.attributes["slice_pfull"]))

            # Add so that it show up on in the title
            sim_var.attributes["short_name"] = "mean $(ClimaAnalysis.short_name(sim_var))"

            # Plot bias
            ClimaAnalysis.Visualize.plot_bias_on_globe!(
                layout,
                sim_var,
                obs_var,
                cmap_extrema = compare_vars_biases_plot_extrema_pfull[short_name],
            )
        end
    end

    # Plot label for the pressure levels
    for (col_idx, p_lvl) in enumerate(real_p_lvls[begin:length(target_p_lvls)])
        CairoMakie.Label(fig_bias_pfull_vars[0, col_idx], "$p_lvl $p_units", tellwidth = false, fontsize = 30)
    end

    # Define figure with at most four columns and the necessary number of rows for
    # all the variables
    num_vars = length(compare_vars_biases_plot_extrema_pfull)
    num_cols = min(4, num_vars)
    num_rows = ceil(Int, num_vars / 4)
    fig_lat_pres = CairoMakie.Figure(size = (650 * num_cols, 450 * num_rows))

    # Initialize ranges for colorbar
    compare_vars_biases_heatmap_extrema_pfull = get_compare_vars_biases_heatmap_extrema_pfull()

    # Take zonal mean and plot lat - pressure heatmap
    for (idx, short_name) in enumerate(keys(sim_obs_comparsion_dict))
        # Compute where the figure should be placed
        # Goes from 1 -> (1,1), 2 -> (1,2), ..., 4 -> (1,4), 5 -> (2,1)
        # for idx -> (col_idx, row_idx)
        row_idx = div(idx - 1, 4) + 1
        col_idx = 1 + rem(idx - 1, 4)
        layout = fig_lat_pres[row_idx, col_idx] = CairoMakie.GridLayout()

        # Load data
        sim_var = sim_obs_comparsion_dict[short_name].sim
        obs_var = sim_obs_comparsion_dict[short_name].obs

        # Take zonal mean (average along lon arithmetically)
        sim_var = ClimaAnalysis.average_lon(sim_var)
        obs_var = ClimaAnalysis.average_lon(obs_var)

        # Compute bias and set short name, long name, and units
        bias_var = sim_var - obs_var
        bias_var = ClimaAnalysis.set_units(bias_var, ClimaAnalysis.units(sim_var))
        bias_var.attributes["short_name"] = "sim-obs_$(ClimaAnalysis.short_name(sim_var))"
        bias_var.attributes["long_name"] = "SIM-OBS_$(ClimaAnalysis.long_name(sim_var))"
        ClimaAnalysis.Visualize.heatmap2D!(
            layout,
            bias_var,
            more_kwargs = Dict(
                :axis => Dict([:yreversed => true]),
                :plot => Dict(
                    :colormap => :vik,
                    :colorrange => compare_vars_biases_heatmap_extrema_pfull[short_name],
                    :colormap => CairoMakie.cgrad(:vik, 11, categorical = true),
                ),
            ),
        )
    end

    # Save figures
    CairoMakie.save(joinpath(leaderboard_base_path, "bias_vars_in_pfull.png"), fig_bias_pfull_vars)
    CairoMakie.save(joinpath(leaderboard_base_path, "bias_lat_pfull_heatmaps.png"), fig_lat_pres)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        error("Usage: julia leaderboard.jl <Filepath to save leaderboard and bias plots> <Filepath to simulation data>")
    end
    leaderboard_base_path = ARGS[begin]
    diagnostics_folder_path = ARGS[begin + 1]
    compute_leaderboard(leaderboard_base_path, diagnostics_folder_path, 3)
    compute_pfull_leaderboard(leaderboard_base_path, diagnostics_folder_path, 6)
end
