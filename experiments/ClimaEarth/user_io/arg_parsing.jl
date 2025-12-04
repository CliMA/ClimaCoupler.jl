import YAML

mode_name_dict = Dict(
    "amip" => AMIPMode,
    "cmip" => CMIPMode,
    "slabplanet" => SlabplanetMode,
    "slabplanet_aqua" => SlabplanetAquaMode,
    "slabplanet_terra" => SlabplanetTerraMode,
    "subseasonal" => SubseasonalMode,
)

"""
    get_coupler_config_dict(config_file)

Read in the configuration file and job ID from the command line.
A dictionary is constructed from the input configuration file and returned.

Since the atmosphere model also uses a configuration file, we read in the atmosphere
configuration file specified in the coupler configuration file (if any), and overwrite
it with the coupler configuration.

The order of priority for overwriting configuration options from lowest to highest is:
    1. ClimaAtmos defaults
    2. ClimaCoupler defaults (defined in `experiments/ClimaEarth/cli_options.jl`)
    3. Command line arguments provided to ClimaCoupler
    4. ClimaAtmos configuration file (if specified in coupler config file)
    5. ClimaCoupler configuration file

# Returns
- `config_dict`: A dictionary mapping configuration keys to the specified settings
"""
function get_coupler_config_dict(config_file)
    # Get the coupler default configuration dictionary, overwritten by any command line arguments
    # Typically the command line arguments are only `config_file` and `job_id`
    coupler_default_cli = parse_commandline(argparse_settings())

    # Extract the job ID from the command line arguments, or from the config file name if not provided
    job_id = coupler_default_cli["job_id"]
    coupler_default_cli["job_id"] =
        isnothing(job_id) ? string(split(split(config_file, '/')[end], '.')[1]) : job_id

    # Load the coupler config file into a dictionary
    coupler_config_dict = YAML.load_file(config_file)

    # Get ClimaAtmos default configuration dictionary
    atmos_default = CA.default_config_dict()
    atmos_config_file = merge(coupler_default_cli, coupler_config_dict)["atmos_config_file"]
    if isnothing(atmos_config_file)
        @info "Using Atmos default configuration"

        # Merge the atmos default config, coupler default config + command line inputs,
        #  and user-provided coupler config
        config_dict = merge(atmos_default, coupler_default_cli, coupler_config_dict)
    else
        @info "Using Atmos configuration from ClimaCoupler in $atmos_config_file"
        atmos_config_dict =
            YAML.load_file(joinpath(pkgdir(ClimaCoupler), atmos_config_file))

        # Merge the atmos default config, coupler default config + command line inputs,
        #  user-provided atmos config, and user-provided coupler config
        config_dict = merge(
            atmos_default,
            coupler_default_cli,
            atmos_config_dict,
            coupler_config_dict,
        )
    end

    # Select the correct timestep for each component model based on which are available
    parse_component_dts!(config_dict)

    return config_dict
end

"""
    get_coupler_args(config_dict)

Extract the necessary arguments from the coupled configuration dictionary.
This function may modify the input dictionary to remove unnecessary keys.

# Arguments
- `config_dict`: A dictionary mapping configuration keys to the specified settings

# Returns
- All arguments needed for the coupled simulation
"""
function get_coupler_args(config_dict::Dict)
    # Vector of TOML files containing model parameters
    # We need to modify this Dict entry to be consistent with ClimaAtmos TOML files
    config_dict["coupler_toml"] = map(config_dict["coupler_toml"]) do file
        isfile(file) ? file : joinpath(pkgdir(ClimaCoupler), file)
    end
    parameter_files = config_dict["coupler_toml"]

    # Make a copy so that we don't modify the original input
    config_dict = copy(config_dict)

    # Simulation-identifying information; Print `config_dict` if requested
    config_dict["print_config_dict"] && @info(config_dict)
    job_id = config_dict["job_id"]
    mode_name = config_dict["mode_name"]
    sim_mode = mode_name_dict[mode_name]
    use_itime = config_dict["use_itime"]

    # Computational simulation setup information
    random_seed = config_dict["unique_seed"] ? time_ns() : 1234
    FT = config_dict["FLOAT_TYPE"] == "Float64" ? Float64 : Float32

    # Time information
    t_end = Float64(Utilities.time_to_seconds(config_dict["t_end"]))
    t_start = Float64(Utilities.time_to_seconds(config_dict["t_start"]))
    start_date = Dates.DateTime(config_dict["start_date"], Dates.dateformat"yyyymmdd")
    Δt_cpl = Float64(Utilities.time_to_seconds(config_dict["dt_cpl"]))

    if use_itime
        t_end = ITime(t_end, epoch = start_date)
        t_start = ITime(t_start, epoch = start_date)
        Δt_cpl = ITime(Δt_cpl, epoch = start_date)
        times = promote(
            t_end,
            t_start,
            Δt_cpl,
            ITime.(values(config_dict["component_dt_dict"]))...,
        )
        t_end, t_start, Δt_cpl = (times[1], times[2], times[3])
        component_dt_dict = Dict(
            component => first(promote(ITime(dt), t_end)) for
            (component, dt) in config_dict["component_dt_dict"]
        )
    else
        component_dt_dict = config_dict["component_dt_dict"]
    end
    # Save solution to integrator.sol at the beginning and end
    saveat = [t_start, t_end]

    # Space information
    share_surface_space = config_dict["share_surface_space"]

    # Checkpointing information
    checkpoint_dt = config_dict["checkpoint_dt"]

    # Restart information
    detect_restart_files = config_dict["detect_restart_files"]
    restart_dir = config_dict["restart_dir"]
    restart_t = config_dict["restart_t"]
    restart_cache = config_dict["restart_cache"]
    save_cache = config_dict["save_cache"]

    # Diagnostics information
    use_coupler_diagnostics = config_dict["use_coupler_diagnostics"]
    use_land_diagnostics = config_dict["use_land_diagnostics"]
    (_, diagnostics_dt) = get_diag_period(t_start, t_end)

    # Physical simulation information
    evolving_ocean = config_dict["evolving_ocean"]

    # Conservation information
    energy_check = config_dict["energy_check"]
    conservation_softfail = config_dict["conservation_softfail"]
    rmse_check = config_dict["rmse_check"]

    # Output information
    output_dir_root = joinpath(config_dict["coupler_output_dir"], job_id)

    # ClimaLand-specific information
    land_model = config_dict["land_model"]
    land_temperature_anomaly = lowercase(config_dict["land_temperature_anomaly"])
    use_land_diagnostics = config_dict["use_land_diagnostics"]
    land_spun_up_ic = config_dict["land_spun_up_ic"]
    bucket_albedo_type = config_dict["bucket_albedo_type"]
    bucket_initial_condition = config_dict["bucket_initial_condition"]

    # Initial condition setting
    era5_initial_condition_dir = config_dict["era5_initial_condition_dir"]

    # Ice model-specific information
    ice_model = config_dict["ice_model"]

    return (;
        job_id,
        sim_mode,
        random_seed,
        FT,
        t_end,
        t_start,
        start_date,
        Δt_cpl,
        component_dt_dict,
        share_surface_space,
        saveat,
        checkpoint_dt,
        detect_restart_files,
        restart_dir,
        restart_t,
        restart_cache,
        save_cache,
        use_coupler_diagnostics,
        diagnostics_dt,
        evolving_ocean,
        energy_check,
        conservation_softfail,
        rmse_check,
        output_dir_root,
        land_model,
        land_temperature_anomaly,
        land_spun_up_ic,
        use_land_diagnostics,
        bucket_albedo_type,
        bucket_initial_condition,
        parameter_files,
        era5_initial_condition_dir,
        ice_model,
    )
end

### Helper functions used in argument parsing ###

"""
    get_diag_period()

Determine the frequency at which to average and output diagnostics based on the
simulation start and end times.

The default periods are:
- 1 month for simulations longer than 90 days
- 10 days for simulations longer than 30 days
- 1 day for simulations longer than 1 day
- 1 hour for simulations shorter than 1 day

# Arguments
- `t_start`: The start time of the simulation
- `t_end`: The end time of the simulation

# Returns
- `period`: A String of how often to average and output diagnostics
- `diagnostics_dt`: A DateTime interval representing the period
"""
function get_diag_period(t_start, t_end)
    sim_duration = float(t_end - t_start)
    secs_per_day = 86400
    if sim_duration >= 90 * secs_per_day
        # if duration >= 90 days, take monthly means
        period = "1months"
        diagnostics_dt = Dates.Month(1)
    elseif sim_duration >= 30 * secs_per_day
        # if duration >= 30 days, take means over 10 days
        period = "10days"
        diagnostics_dt = Dates.Day(10)
    elseif sim_duration >= secs_per_day
        # if duration >= 1 day, take daily means
        period = "1days"
        diagnostics_dt = Dates.Day(1)
    else
        # if duration < 1 day, take hourly means
        period = "1hours"
        diagnostics_dt = Dates.Hour(1)
    end
    return (period, diagnostics_dt)
end

"""
    parse_component_dts!(config_dict)

Check which timesteps are specified in the config file, and use them to choose
the correct timestep for each component model.
If all component timesteps `dt_\$component` are specified in the config file, use those
and remove `dt` if it was provided.
Otherwise, use the generic component timestep `dt` specified in the config file.
If some (but not all) component timesteps and the generic timestep `dt` are specified,
use the generic timestep and remove the others from the config dict.

The timestep for each component model is stored in the `component_dt_dict` field of the config dict.

# Arguments
- `config_dict`: A dictionary mapping configuration keys to the specified settings
"""
function parse_component_dts!(config_dict)
    # Retrieve coupling timestep
    Δt_cpl = Float64(Utilities.time_to_seconds(config_dict["dt_cpl"]))

    # Specify component model names
    component_dt_names = ["dt_atmos", "dt_land", "dt_ocean", "dt_seaice"]
    component_dt_dict = Dict{String, typeof(Δt_cpl)}()
    # check if all component dt's are specified
    if all(key -> !isnothing(config_dict[key]), component_dt_names)
        # when all component dt's are specified, ignore the dt field
        if haskey(config_dict, "dt")
            @warn "Removing dt in favor of individual component dt's"
            delete!(config_dict, "dt")
        end
        for key in component_dt_names
            component_dt = Float64(Utilities.time_to_seconds(config_dict[key]))
            @assert isapprox(Δt_cpl % component_dt, 0.0) "Coupler dt must be divisible by all component dt's\n dt_cpl = $Δt_cpl\n $key = $component_dt"
            component_dt_dict[key] = component_dt
        end
    else
        # when not all component dt's are specified, use the dt field
        @assert haskey(config_dict, "dt") "dt or (dt_atmos, dt_land, dt_ocean, and dt_seaice) must be specified"
        for key in component_dt_names
            if !isnothing(config_dict[key])
                @warn "Removing $key from config in favor of dt because not all component dt's are specified"
            end
            delete!(config_dict, key)
            component_dt_dict[key] = Float64(Utilities.time_to_seconds(config_dict["dt"]))
        end
    end
    config_dict["component_dt_dict"] = component_dt_dict
    return nothing
end
