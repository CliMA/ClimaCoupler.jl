# For generating plots
import ClimaAnalysis
import GeoMakie
import CairoMakie

@info "Error against observations"

# Tuple of short names for loading simulation and observational data
sim_obs_short_names_no_pr = [
    ("rsdt", "solar_mon"),
    ("rsut", "toa_sw_all_mon"),
    ("rlut", "toa_lw_all_mon"),
    ("rsutcs", "toa_sw_clr_t_mon"),
    ("rlutcs", "toa_lw_clr_t_mon"),
    ("rsds", "sfc_sw_down_all_mon"),
    ("rsus", "sfc_sw_up_all_mon"),
    ("rlds", "sfc_lw_down_all_mon"),
    ("rlus", "sfc_lw_up_all_mon"),
    ("rsdscs", "sfc_sw_down_clr_t_mon"),
    ("rsuscs", "sfc_sw_up_clr_t_mon"),
    ("rldscs", "sfc_lw_down_clr_t_mon"),
]

compare_vars_biases_plot_extrema = Dict(
    "pr" => (-5.0, 5.0),
    "rsdt" => (-2.0, 2.0),
    "rsut" => (-50.0, 50.0),
    "rlut" => (-50.0, 50.0),
    "rsutcs" => (-20.0, 20.0),
    "rlutcs" => (-20.0, 20.0),
    "rsds" => (-50.0, 50.0),
    "rsus" => (-10.0, 10.0),
    "rlds" => (-50.0, 50.0),
    "rlus" => (-50.0, 50.0),
    "rsdscs" => (-10.0, 10.0),
    "rsuscs" => (-10.0, 10.0),
    "rldscs" => (-20.0, 20.0),
)


if length(ARGS) < 1
    error("Usage: leaderboard.jl <path of folder containing NetCDF files>")
end

# Path to saved leaderboards
leaderboard_base_path = ARGS[begin]

# Path to simulation data
diagnostics_folder_path = ARGS[begin]

# Dict for loading in simulation data
sim_var_dict = Dict{String,Any}(
    "pr" = () -> begin
        sim_var = get(ClimaAnalysis.SimDir(diagnostics_folder_path), short_name = "pr")
        sim_var = ClimaAnalysis.convert_units(
            sim_var,
            "mm/day",
            conversion_function = x -> x .* Float32(-86400),
        )
        sim_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
        return sim_var
    end,
)

# Loop to load the rest of the simulation data
for (short_name, _) in sim_obs_short_names_no_pr
    sim_var_dict[short_name] =
        () -> begin
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = short_name,
            )
            sim_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end
end

# Dict for loading observational data
obs_var_dict = Dict{String,Any}(
    "pr" =>
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                joinpath(
                    @clima_artifact("precipitation_obs"),
                    "gpcp.precip.mon.mean.197901-202305.nc",
                ),
                "precip",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            return obs_var
        end,
)

# Loop to load the rest of the observational data
for (sim_name, obs_name) in sim_obs_short_names_no_pr
    obs_var_dict[sim_name] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                joinpath(
                    @clima_artifact("radiation_obs"),
                    "CERES_EBAF_Ed4.2_Subset_200003-201910.nc",
                ),
                obs_name,
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            # Convert from W m-2 to W m^-2
            ClimaAnalysis.units(obs_var) == "W m-2" ?
            obs_var = ClimaAnalysis.set_units(obs_var, "W m^-2") :
            error("Did not expect $(ClimaAnalysis.units(obs_var)) for the units")
            return obs_var
        end
end

# Set up dict for storing simulation and observational data after processing
sim_obs_comparsion_dict = Dict()
seasons = ["ANN", "MAM", "JJA", "SON", "DJF"]

# Print dates for debugging
pr_var = sim_var_dict["pr"]() # it shouldn't matter what short name we use
output_dates =
    Dates.DateTime(pr_var.attributes["start_date"]) .+
    Dates.Second.(ClimaAnalysis.times(pr_var))
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
    obs_var_seasons = [obs_var, obs_var_seasons...]
    sim_var_seasons = [sim_var, sim_var_seasons...]

    # Take time average
    obs_var_seasons = obs_var_seasons .|> ClimaAnalysis.average_time
    sim_var_seasons = sim_var_seasons .|> ClimaAnalysis.average_time

    # Add "mean " for plotting the title
    for sim_var in sim_var_seasons
        sim_var.attributes["short_name"] = "mean $(ClimaAnalysis.short_name(sim_var))"
    end

    # Save observation and simulation data
    sim_obs_comparsion_dict[short_name] = Dict(
        season => (sim_var_s, obs_var_s) for
        (season, sim_var_s, obs_var_s) in zip(seasons, sim_var_seasons, obs_var_seasons)
    )
end

compare_vars_biases_groups = [
    ["pr", "rsdt", "rsut", "rlut"],
    ["rsds", "rsus", "rlds", "rlus"],
    ["rsutcs", "rlutcs", "rsdscs", "rsuscs", "rldscs"],
]

# Plot bias plots
for season in seasons
    for compare_vars_biases in compare_vars_biases_groups
        fig_bias = CairoMakie.Figure(; size = (600, 300 * length(compare_vars_biases)))
        for (loc, short_name) in enumerate(compare_vars_biases)
            ClimaAnalysis.Visualize.plot_bias_on_globe!(
                fig_bias,
                sim_obs_comparsion_dict[short_name][season]...,
                cmap_extrema = compare_vars_biases_plot_extrema[short_name],
                p_loc = (loc, 1),
            )
        end
        # Do if and else statement for naming files appropriately
        if season != "ANN"
            CairoMakie.save(
                joinpath(
                    leaderboard_base_path,
                    "bias_$(first(compare_vars_biases))_$season.png",
                ),
                fig_bias,
            )
        else
            CairoMakie.save(
                joinpath(
                    leaderboard_base_path,
                    "bias_$(first(compare_vars_biases))_total.png",
                ),
                fig_bias,
            )
        end
    end
end

# Plot leaderboard
# Load data into RMSEVariables
rmse_var_pr = ClimaAnalysis.read_rmses(
    joinpath(@clima_artifact("cmip_model_rmse"), "pr_rmse_amip_pr_amip_5yr.csv"),
    "pr",
    units = "mm / day",
)
rmse_var_rsut = ClimaAnalysis.read_rmses(
    joinpath(@clima_artifact("cmip_model_rmse"), "rsut_rmse_amip_rsut_amip_5yr.csv"),
    "rsut",
    units = "W m^-2",
)
rmse_var_rlut = ClimaAnalysis.read_rmses(
    joinpath(@clima_artifact("cmip_model_rmse"), "rlut_rmse_amip_rlut_amip_5yr.csv"),
    "rlut",
    units = "W m^-2",
)

# Add models and units for CliMA
rmse_var_pr = ClimaAnalysis.add_model(rmse_var_pr, "CliMA")
ClimaAnalysis.add_unit!(rmse_var_pr, "CliMA", "mm / day")

rmse_var_rsut = ClimaAnalysis.add_model(rmse_var_rsut, "CliMA")
ClimaAnalysis.add_unit!(rmse_var_rsut, "CliMA", "W m^-2")

rmse_var_rlut = ClimaAnalysis.add_model(rmse_var_rlut, "CliMA")
ClimaAnalysis.add_unit!(rmse_var_rlut, "CliMA", "W m^-2")

# Add RMSE for the CliMA model and for each season
for season in seasons
    rmse_var_pr["CliMA", season] =
        ClimaAnalysis.global_rmse(sim_obs_comparsion_dict["pr"][season]...)
    rmse_var_rsut["CliMA", season] =
        ClimaAnalysis.global_rmse(sim_obs_comparsion_dict["rsut"][season]...)
    rmse_var_rlut["CliMA", season] =
        ClimaAnalysis.global_rmse(sim_obs_comparsion_dict["rlut"][season]...)
end

# Plot box plots
rmse_vars = (rmse_var_pr, rmse_var_rsut, rmse_var_rlut)
fig_leaderboard = CairoMakie.Figure(; size = (800, 300 * 3 + 400), fontsize = 20)
for (loc, rmse_var) in enumerate(rmse_vars)
    ClimaAnalysis.Visualize.plot_boxplot!(
        fig_leaderboard,
        rmse_var,
        ploc = (loc, 1),
        best_and_worst_category_name = "ANN",
    )
end

# Plot leaderboard
ClimaAnalysis.Visualize.plot_leaderboard!(
    fig_leaderboard,
    rmse_vars...,
    best_category_name = "ANN",
    ploc = (4, 1),
)
CairoMakie.save(joinpath(leaderboard_base_path, "bias_leaderboard.png"), fig_leaderboard)
