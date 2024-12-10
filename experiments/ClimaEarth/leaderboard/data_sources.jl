import ClimaAnalysis
import ClimaUtilities.ClimaArtifacts: @clima_artifact

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

"""
    get_compare_vars_biases_groups()

Return an array of arrays of short names.

This determines which variables are grouped on the bias plots.
"""
function get_compare_vars_biases_groups()
    compare_vars_biases_groups = [
        ["pr", "rsdt", "rsut", "rlut"],
        ["rsds", "rsus", "rlds", "rlus"],
        ["rsutcs", "rlutcs", "rsdscs", "rsuscs", "rldscs"],
    ]
    return compare_vars_biases_groups
end

"""
    get_compare_vars_biases_plot_extrema()

Return a dictionary mapping short names to ranges for the bias plots. This is
used by the function `compute_leaderboard`.

To add a variable to the leaderboard, add a key-value pair to the dictionary
`compare_vars_biases_plot_extrema` whose key is a short name key is the same
short name in `sim_var_dict` and the value is a tuple, where the first element
is the lower bound and the last element is the upper bound for the bias plots.
"""
function get_compare_vars_biases_plot_extrema()
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
    return compare_vars_biases_plot_extrema
end

"""
    get_sim_var_dict(diagnostics_folder_path)

Return a dictionary mapping short names to `OutputVar` containing preprocessed
simulation data. This is used by the function `compute_leaderboard`.

To add a variable for the leaderboard, add a key-value pair to the dictionary
`sim_var_dict` whose key is the short name of the variable and the value is an
anonymous function that returns a `OutputVar`. For each variable, any
preprocessing should be done in the corresponding anonymous function which
includes unit conversion and shifting the dates.

The variable should have only three dimensions: latitude, longitude, and time.
"""
function get_sim_var_dict(diagnostics_folder_path)
    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}(
        "pr" =>
            () -> begin
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

    # Add "pr" and the necessary preprocessing
    for (short_name, _) in sim_obs_short_names_no_pr
        sim_var_dict[short_name] =
            () -> begin
                sim_var = get(ClimaAnalysis.SimDir(diagnostics_folder_path), short_name = short_name)
                sim_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
                return sim_var
            end
    end
    return sim_var_dict
end

"""
    get_obs_var_dict()

Return a dictionary mapping short names to `OutputVar` containing preprocessed
observational data. This is used by the function `compute_leaderboard`.

To add a variable for the leaderboard, add a key-value pair to the dictionary
`obs_var_dict` whose key is the short name of the variable and the value is an
anonymous function that returns a `OutputVar`. The function must take in a
start date which is used to align the times in the observational data to match
the simulation data. The short name must be the same as in `sim_var_dict` in the
function `sim_var_dict`. For each variable, any preprocessing is done in the
corresponding anonymous function which includes unit conversion and shifting the
dates.

The variable should have only three dimensions: latitude, longitude, and time.
"""
function get_obs_var_dict()
    # Add "pr" and the necessary preprocessing
    obs_var_dict = Dict{String, Any}(
        "pr" =>
            (start_date) -> begin
                obs_var = ClimaAnalysis.OutputVar(
                    joinpath(@clima_artifact("precipitation_obs"), "precip.mon.mean.nc"),
                    "precip",
                    new_start_date = start_date,
                    shift_by = Dates.firstdayofmonth,
                )
                return obs_var
            end,
    )

    # Loop to load the rest of the observational data and the necessary preprocessing
    for (sim_name, obs_name) in sim_obs_short_names_no_pr
        obs_var_dict[sim_name] =
            (start_date) -> begin
                obs_var = ClimaAnalysis.OutputVar(
                    joinpath(@clima_artifact("radiation_obs"), "CERES_EBAF_Ed4.2_Subset_200003-201910.nc"),
                    obs_name,
                    new_start_date = start_date,
                    shift_by = Dates.firstdayofmonth,
                )
                # Convert from W m-2 to W m^-2
                ClimaAnalysis.units(obs_var) == "W m-2" ? obs_var = ClimaAnalysis.set_units(obs_var, "W m^-2") :
                error("Did not expect $(ClimaAnalysis.units(obs_var)) for the units")
                return obs_var
            end
    end
    return obs_var_dict
end

"""
    get_rmse_var_dict()

Return a dictionary mapping short names to `RMSEVariable` containing information about
RMSEs of models for the short name of the variable. This is used by the function
`compute_leaderboard`.
"""
function get_rmse_var_dict()
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

    # Store the rmse vars in a dictionary
    rmse_var_dict = ClimaAnalysis.OrderedCollections.OrderedDict(
        "pr" => rmse_var_pr,
        "rsut" => rmse_var_rsut,
        "rlut" => rmse_var_rlut,
    )
    return rmse_var_dict
end

"""
    get_sim_var_in_pfull_dict(diagnostics_folder_path)

Return a dictionary mapping short names to `OutputVar` containing preprocessed
simulation data in pressure space. This is used by `compute_pfull_leaderboard`.

To add a variable for the leaderboard, add a key-value pair to the dictionary
`sim_var_dict` whose key is the short name of the variable and the value is an
anonymous function that returns a `OutputVar`. For each variable, any
preprocessing should be done in the corresponding anonymous function which
includes unit conversion and shifting the dates.

The variable should have only four dimensions: latitude, longitude, time, and
pressure. The units of pressure should be in hPa.
"""
function get_sim_var_in_pfull_dict(diagnostics_folder_path)
    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    sim_var_pfull_dict = Dict{String, Any}()

    short_names = ["ta", "hur", "hus"]
    for short_name in short_names
        sim_var_pfull_dict[short_name] =
            () -> begin
                sim_var = get(simdir, short_name = short_name)
                pfull_var = get(simdir, short_name = "pfull")

                (ClimaAnalysis.units(sim_var) == "") && (sim_var = ClimaAnalysis.set_units(sim_var, "unitless"))
                (ClimaAnalysis.units(sim_var) == "kg kg^-1") &&
                    (sim_var = ClimaAnalysis.set_units(sim_var, "unitless"))

                sim_in_pfull_var = ClimaAnalysis.Atmos.to_pressure_coordinates(sim_var, pfull_var)
                sim_in_pfull_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_in_pfull_var)
                sim_in_pfull_var = ClimaAnalysis.convert_dim_units(
                    sim_in_pfull_var,
                    "pfull",
                    "hPa";
                    conversion_function = x -> 0.01 * x,
                )
                return sim_in_pfull_var
            end
    end
    return sim_var_pfull_dict
end

"""
    get_obs_var_in_pfull_dict()

Return a dictionary mapping short names to `OutputVar` containing preprocessed
observational data in pressure coordinates. This is used by
`compute_pfull_leaderboard`.

To add a variable for the leaderboard, add a key-value pair to the dictionary
`obs_var_dict` whose key is the short name of the variable and the value is an anonymous
function that returns a `OutputVar`. The function must take in a start date
which is used to align the times in the observational data to match the
simulation data. The short name must be the same as in `sim_var_pfull_dict` in
the function `get_sim_var_in_pfull_dict`. Any preprocessing is done in the
function which includes unit conversion and shifting the dates.

The variable should have only four dimensions: latitude, longitude, time, and
pressure. The units of pressure should be in hPa.
"""
function get_obs_var_in_pfull_dict()
    artifact_path = joinpath(
        @clima_artifact("era5_monthly_averages_pressure_levels_1979_2024"),
        "era5_monthly_averages_pressure_levels_197901-202410.nc",
    )
    short_names_pairs = [("ta", "t"), ("hus", "q")]
    obs_var_dict = Dict{String, Any}()
    for (short_name, era5_short_name) in short_names_pairs
        obs_var_dict[short_name] =
            (start_date) -> begin
                obs_var = ClimaAnalysis.OutputVar(
                    artifact_path,
                    era5_short_name,
                    new_start_date = start_date,
                    shift_by = Dates.firstdayofmonth,
                )
                (ClimaAnalysis.units(obs_var) == "kg kg**-1") &&
                    (obs_var = ClimaAnalysis.set_units(obs_var, "unitless"))
                obs_var = ClimaAnalysis.Var.convert_dim_units(
                    obs_var,
                    "pressure_level",
                    "hPa";
                    conversion_function = x -> 0.01 * x,
                )
                return obs_var
            end
    end

    obs_var_dict["hur"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                artifact_path,
                "r",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            obs_var = ClimaAnalysis.Var.convert_dim_units(
                obs_var,
                "pressure_level",
                "hPa";
                conversion_function = x -> 0.01 * x,
            )
            # Convert from percentages (e.g. 120%) to decimal (1.20)
            obs_var = ClimaAnalysis.Var.convert_units(obs_var, "unitless", conversion_function = x -> 0.01 * x)
            return obs_var
        end
    return obs_var_dict
end

"""
    get_compare_vars_biases_plot_extrema_pfull()

Return a dictionary mapping short names to ranges for the bias plots. This is
used by the function `compute_pfull_leaderboard`.

To add a variable to the leaderboard, add a key-value pair to the dictionary
`compare_vars_biases_plot_extrema` whose key is a short name key is the same
short name in `sim_var_pfull_dict` in the function `get_sim_var_pfull_dict` and
the value is a tuple, where the first element is the lower bound and the last
element is the upper bound for the bias plots.
"""
function get_compare_vars_biases_plot_extrema_pfull()
    compare_vars_biases_plot_extrema = Dict("ta" => (-5.0, 5.0), "hur" => (-0.4, 0.4), "hus" => (-0.0015, 0.0015))
    return compare_vars_biases_plot_extrema
end

"""
    get_compare_vars_biases_heatmap_extrema_pfull()

Return a dictionary mapping short names to ranges for the heatmaps. This is
used by the function `compute_pfull_leaderboard`.

To add a variable to the leaderboard, add a key-value pair to the dictionary
`compare_vars_biases_heatmap_extrema` whose key is a short name key is the same
short name in `sim_var_pfull_dict` in the function `get_sim_var_pfull_dict` and
the value is a tuple, where the first element is the lower bound and the last
element is the upper bound for the bias plots.
"""
function get_compare_vars_biases_heatmap_extrema_pfull()
    compare_vars_biases_heatmap_extrema = Dict("ta" => (-10.0, 10.0), "hur" => (-0.4, 0.4), "hus" => (-0.001, 0.001))
    return compare_vars_biases_heatmap_extrema
end
