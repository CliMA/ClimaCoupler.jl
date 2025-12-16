import ClimaAnalysis
import ClimaUtilities.ClimaArtifacts: @clima_artifact

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
        ["hfss", "hfls"],
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
        "hfss" => (-50.0, 50.0),
        "hfls" => (-50.0, 50.0),
    )
    return compare_vars_biases_plot_extrema
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
    rmse_var_dict = ClimaAnalysis.Var.OrderedDict(
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

Some configurations have pressure inversions at low elevations, which prevents the
conversion to pressure coordiantes. This function discards everything below 80 meters of
elevation to avoid that.
"""
function get_sim_var_in_pfull_dict(diagnostics_folder_path)
    available_short_names =
        ClimaAnalysis.available_vars(ClimaAnalysis.SimDir(diagnostics_folder_path))
    sim_var_pfull_dict = Dict{String, Any}()

    short_names = get_short_names_monthly_averages(diagnostics_folder_path)
    available_short_names = intersect(short_names, Set(["ta", "hur", "hus"]))
    for short_name in short_names
        short_name in available_short_names && (
            sim_var_pfull_dict[short_name] =
                () -> begin
                    sim_var = get(
                        ClimaAnalysis.SimDir(diagnostics_folder_path),
                        short_name = short_name,
                        reduction = "average",
                        period = "1M",
                    )
                    pfull_var = get(
                        ClimaAnalysis.SimDir(diagnostics_folder_path),
                        short_name = "pfull",
                        reduction = "average",
                        period = "1M",
                    )

                    (ClimaAnalysis.units(sim_var) == "") &&
                        (sim_var = ClimaAnalysis.set_units(sim_var, "unitless"))
                    (ClimaAnalysis.units(sim_var) == "kg kg^-1") &&
                        (sim_var = ClimaAnalysis.set_units(sim_var, "unitless"))

                    # For certain grid configurations, certain columns can have pressure
                    # inversions at low elevation, which prevents the conversion between
                    # pressure and altitude. To avoid that, we exclude everything below 80m.
                    # NOTE: This was empirically found.
                    pfull_var_windowed =
                        ClimaAnalysis.window(pfull_var, "z", left = 80)
                    sim_var_windowed = ClimaAnalysis.window(sim_var, "z", left = 80)

                    sim_in_pfull_var = ClimaAnalysis.Atmos.to_pressure_coordinates(
                        sim_var_windowed,
                        pfull_var_windowed,
                    )
                    sim_in_pfull_var =
                        ClimaAnalysis.shift_to_start_of_previous_month(sim_in_pfull_var)
                    sim_in_pfull_var = ClimaAnalysis.convert_dim_units(
                        sim_in_pfull_var,
                        "pfull",
                        "hPa";
                        conversion_function = x -> 0.01 * x,
                    )
                    return sim_in_pfull_var
                end
        )
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
            obs_var = ClimaAnalysis.Var.convert_units(
                obs_var,
                "unitless",
                conversion_function = x -> 0.01 * x,
            )
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
    compare_vars_biases_plot_extrema =
        Dict("ta" => (-5.0, 5.0), "hur" => (-0.4, 0.4), "hus" => (-0.0015, 0.0015))
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
    compare_vars_biases_heatmap_extrema =
        Dict("ta" => (-10.0, 10.0), "hur" => (-0.4, 0.4), "hus" => (-0.001, 0.001))
    return compare_vars_biases_heatmap_extrema
end

"""
    get_short_names_of_monthly_averages(diagnostics_folder_path)

Get all the short names of the monthly averages.
"""
function get_short_names_monthly_averages(diagnostics_folder_path)
    available_short_names = Set{String}()
    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in ClimaAnalysis.available_vars(simdir)
        for reduction in ClimaAnalysis.available_reductions(simdir, short_name = short_name)
            for period in ClimaAnalysis.available_periods(
                simdir,
                short_name = short_name,
                reduction = reduction,
            )
                if reduction == "average" && period == "1M"
                    push!(available_short_names, short_name)
                end
            end
        end
    end
    return available_short_names
end
