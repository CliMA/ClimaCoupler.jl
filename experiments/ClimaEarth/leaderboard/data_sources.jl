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

Return a dictionary mapping short names to ranges for the bias plots.

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
simulation data.

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
observational data.

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
RMSEs of models for the short name of the variable.
"""
function get_rmse_var_dict()
    # Load data into RMSEVariables
    rmse_var_names = ["pr", "rsut", "rlut"]
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
    rmse_var_dict = Dict("pr" => rmse_var_pr, "rsut" => rmse_var_rsut, "rlut" => rmse_var_rlut)
    return rmse_var_dict
end
