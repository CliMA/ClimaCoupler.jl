import ClimaAnalysis
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Dates

export get_sim_var_dict, get_obs_var_dict

# Tuple of short names for loading simulation and observational data
const sim_obs_short_names_rad = [
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

const sim_obs_short_names_sf = [("hfss", "msshf"), ("hfls", "mslhf")]

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

"""
    get_sim_var_dict(diagnostics_folder_path)

Return a dictionary mapping short names to `OutputVar` containing preprocessed
simulation data. This is used by the function `compute_leaderboard`.

To add a variable for the leaderboard, add a key-value pair to the dictionary
`sim_var_dict` whose key is the short name of the variable and the value is an
anonymous function that returns a `OutputVar`. For each variable, any
preprocessing should be done in the corresponding anonymous function which
includes unit conversion.

The variable should have only three dimensions: latitude, longitude, and time.
"""
function get_sim_var_dict(diagnostics_folder_path)
    available_short_names = get_short_names_monthly_averages(diagnostics_folder_path)
    # Dict for loading in simulation data
    sim_var_dict = Dict{String, Any}()
    # Add "pr" and the necessary preprocessing
    "pr" in available_short_names && (
        sim_var_dict["pr"] =
            () -> begin
                sim_var = get(
                    ClimaAnalysis.SimDir(diagnostics_folder_path),
                    short_name = "pr",
                    reduction = "average",
                    period = "1M",
                )
                sim_var = ClimaAnalysis.convert_units(
                    sim_var,
                    "mm/day",
                    conversion_function = x -> x .* Float32(-86400),
                )
                # For ClimaDiagnostics v0.3 and later, the dates are saved at
                # the start of the reduction period
                pkgversion(CD) < v"0.3" && (
                    sim_var = ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
                )
                return sim_var
            end
    )

    for (short_name, _) in vcat(sim_obs_short_names_rad, sim_obs_short_names_sf)
        short_name in available_short_names && (
            sim_var_dict[short_name] =
                () -> begin
                    sim_var = get(
                        ClimaAnalysis.SimDir(diagnostics_folder_path),
                        short_name = short_name,
                        reduction = "average",
                        period = "1M",
                    )
                    # For ClimaDiagnostics v0.3 and later, the dates are saved at
                    # the start of the reduction period
                    pkgversion(CD) < v"0.3" && (
                        sim_var =
                            ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
                    )
                    return sim_var
                end
        )
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
    obs_var_dict = Dict{String, Any}()
    obs_var_dict["pr"] =
        (start_date) -> begin
            obs_var = ClimaAnalysis.OutputVar(
                joinpath(@clima_artifact("precipitation_obs"), "precip.mon.mean.nc"),
                "precip",
                new_start_date = start_date,
                shift_by = Dates.firstdayofmonth,
            )
            return obs_var
        end

    for (sim_name, obs_name) in sim_obs_short_names_sf
        obs_var_dict[sim_name] =
            (start_date) -> begin
                obs_var = ClimaAnalysis.OutputVar(
                    joinpath(
                        @clima_artifact(
                            "era5_monthly_averages_surface_single_level_1979_2024"
                        ),
                        "era5_monthly_averages_surface_single_level_197901-202410.nc",
                    ),
                    obs_name,
                    new_start_date = start_date,
                    shift_by = Dates.firstdayofmonth,
                )
                (ClimaAnalysis.units(obs_var) == "W m**-2") && (
                    obs_var = ClimaAnalysis.convert_units(
                        obs_var,
                        "W m^-2",
                        conversion_function = units -> units * -1.0,
                    )
                )
                return obs_var
            end
    end

    for (sim_name, obs_name) in sim_obs_short_names_rad
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
    return obs_var_dict
end
