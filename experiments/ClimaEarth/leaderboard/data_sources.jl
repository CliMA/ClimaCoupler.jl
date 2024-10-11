import ClimaAnalysis
import ClimaUtilities.ClimaArtifacts: @clima_artifact

# For plotting bias plots
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

# To add a variable to the bias plots, add the variable to sim_var_dict, obs_var_dict,
# compare_vars_biases_groups, and compare_vars_biases_plot_extrema
# To add a variable to the leaderboard, add the variable to rmse_var_names and rmse_var_dict

# Dict for loading in simulation data
sim_var_dict = Dict{String, Any}(
    "pr" =>
        () -> begin
            sim_var = get(ClimaAnalysis.SimDir(diagnostics_folder_path), short_name = "pr")
            sim_var =
                ClimaAnalysis.convert_units(sim_var, "mm/day", conversion_function = x -> x .* Float32(-86400))
            sim_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end,
)

# Loop to load the rest of the simulation data
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

# Add "pr" and the necessary preprocessing
for (short_name, _) in sim_obs_short_names_no_pr
    sim_var_dict[short_name] =
        () -> begin
            sim_var = get(ClimaAnalysis.SimDir(diagnostics_folder_path), short_name = short_name)
            sim_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end
end

# Dict for loading observational data
# Add "pr" and the necessary preprocessing
obs_var_dict = Dict{String, Any}(
    "pr" =>
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                joinpath(@clima_artifact("precipitation_obs"), "gpcp.precip.mon.mean.197901-202305.nc"),
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

# Used to organize plots
compare_vars_biases_groups = [
    ["pr", "rsdt", "rsut", "rlut"],
    ["rsds", "rsus", "rlds", "rlus"],
    ["rsutcs", "rlutcs", "rsdscs", "rsuscs", "rldscs"],
]

# For plotting box plots and leaderboard

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
