import YAML

mode_name_dict = Dict(
    "amip" => AMIPMode,
    "slabplanet" => SlabplanetMode,
    "slabplanet_aqua" => SlabplanetAquaMode,
    "slabplanet_terra" => SlabplanetTerraMode,
    "slabplanet_eisenman" => SlabplanetEisenmanMode,
)

"""
    get_coupler_config()

Read in the configuration file and job ID from the command line.
A dictionary is constructed from the input configuration file and returned.

# Returns
- `config_dict`: A dictionary mapping configuration keys to the specified settings
"""
function get_coupler_config()
    # Read in command line arguments
    parsed_args = parse_commandline(argparse_settings())

    # Extract the configuration file and job ID
    config_file = parsed_args["config_file"]
    job_id = parsed_args["job_id"]
    # Get the job ID from the config file string if not provided
    job_id = isnothing(job_id) ? string(split(split(config_file, '/')[end], '.')[1]) : job_id

    # Read in config dictionary from file, overriding the defaults in `parsed_args`
    config_dict = merge(parsed_args, YAML.load_file(parsed_args["config_file"]))
    config_dict["job_id"] = job_id
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
    # Simulation-identifying information; Print `config_dict` if requested
    config_dict["print_config_dict"] && @info(config_dict)
    job_id = config_dict["job_id"]
    mode_name = config_dict["mode_name"]
    sim_mode = mode_name_dict[mode_name]

    # Computational simulation setup information
    random_seed = config_dict["unique_seed"] ? time_ns() : 1234
    FT = config_dict["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
    comms_ctx = Utilities.get_comms_context(config_dict)

    # Time information
    t_end = Float64(Utilities.time_to_seconds(config_dict["t_end"]))
    t_start = Float64(Utilities.time_to_seconds(config_dict["t_start"]))
    date0 = date = Dates.DateTime(config_dict["start_date"], Dates.dateformat"yyyymmdd")
    Δt_cpl = Float64(Utilities.time_to_seconds(config_dict["dt_cpl"]))
    saveat = Float64(Utilities.time_to_seconds(config_dict["dt_save_to_sol"]))
    component_dt_dict = config_dict["component_dt_dict"]

    # Checkpointing information
    hourly_checkpoint = config_dict["hourly_checkpoint"]
    hourly_checkpoint_dt = config_dict["hourly_checkpoint_dt"]

    # Restart information
    restart_dir = config_dict["restart_dir"]
    restart_t = Int(config_dict["restart_t"])

    # Diagnostics information
    use_coupler_diagnostics = config_dict["use_coupler_diagnostics"]
    use_land_diagnostics = config_dict["use_land_diagnostics"]
    calendar_dt = config_dict["calendar_dt"]

    # Physical simulation information
    evolving_ocean = config_dict["evolving_ocean"]
    mono_surface = config_dict["mono_surface"]
    turb_flux_partition = config_dict["turb_flux_partition"]

    # Conservation information
    energy_check = config_dict["energy_check"]
    conservation_softfail = config_dict["conservation_softfail"]

    # Output information
    output_dir_root = config_dict["coupler_output_dir"]
    plot_diagnostics = config_dict["plot_diagnostics"]

    # ClimaLand-specific information
    land_domain_type = config_dict["land_domain_type"]
    land_albedo_type = config_dict["land_albedo_type"]
    land_initial_condition = config_dict["land_initial_condition"]
    land_temperature_anomaly = config_dict["land_temperature_anomaly"]
    use_land_diagnostics = config_dict["use_land_diagnostics"]

    return (;
        job_id,
        sim_mode,
        random_seed,
        FT,
        comms_ctx,
        t_end,
        t_start,
        date0,
        date,
        Δt_cpl,
        component_dt_dict,
        saveat,
        hourly_checkpoint,
        hourly_checkpoint_dt,
        restart_dir,
        restart_t,
        use_coupler_diagnostics,
        calendar_dt,
        evolving_ocean,
        mono_surface,
        turb_flux_partition,
        energy_check,
        conservation_softfail,
        output_dir_root,
        plot_diagnostics,
        land_domain_type,
        land_albedo_type,
        land_initial_condition,
        land_temperature_anomaly,
        use_land_diagnostics,
    )
end

"""
    get_atmos_args(atmos_config_dict)

Extract the necessary arguments from the atmosphere configuration dictionary.

# Arguments
- `atmos_config_dict`: A dictionary mapping atmosphere configuration keys to the specified settings

# Returns
- All arguments needed for the atmosphere simulation
"""
function get_atmos_args(atmos_config_dict)
    dt_rad = atmos_config_dict["dt_rad"]
    output_default_diagnostics = atmos_config_dict["output_default_diagnostics"]

    return (; dt_rad, output_default_diagnostics)
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
- `calendar_dt`: A DateTime interval representing the period
"""
function get_diag_period(t_start, t_end)
    sim_duration = t_end - t_start
    secs_per_day = 86400
    if sim_duration >= 90 * secs_per_day
        # if duration >= 90 days, take monthly means
        period = "1months"
        calendar_dt = Dates.Month(1)
    elseif sim_duration >= 30 * secs_per_day
        # if duration >= 30 days, take means over 10 days
        period = "10days"
        calendar_dt = Dates.Day(10)
    elseif sim_duration >= secs_per_day
        # if duration >= 1 day, take daily means
        period = "1days"
        calendar_dt = Dates.Day(1)
    else
        # if duration < 1 day, take hourly means
        period = "1hours"
        calendar_dt = Dates.Hour(1)
    end
    return (period, calendar_dt)
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
    component_dt_dict = Dict{String, Float64}()
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

"""
    add_extra_diagnostics!(config_dict)

Conditionally add extra diagnostics to the config dictionary based on the
simulation type and flag to use diagnostics. Currently, the only extra
diagnostic is the atmosphere TOA net flux for AMIP simulations, but more diagnostics
can be added for any component by following the structure in `climaatmos_extra_diags.jl`.

The added atmosphere diagnostics are added to the `extra_atmos_diagnostics` field
of the config dict. The `calendar_dt` field is also added to the config dict to
coordinate the output frequency of the diagnostics.

# Arguments
- `config_dict`: A dictionary mapping configuration keys to the specified settings
"""
function add_extra_diagnostics!(config_dict)
    # Diagnostics information
    mode_name = config_dict["mode_name"]
    use_coupler_diagnostics = config_dict["use_coupler_diagnostics"]
    t_end = Float64(Utilities.time_to_seconds(config_dict["t_end"]))
    t_start = Float64(Utilities.time_to_seconds(config_dict["t_start"]))
    calendar_dt = nothing
    if mode_name == "amip" && use_coupler_diagnostics
        @info "Using default AMIP diagnostics"
        (period, calendar_dt) = get_diag_period(t_start, t_end)

        # Additional atmosphere diagnostics
        !haskey(config_dict, "extra_atmos_diagnostics") &&
            (config_dict["extra_atmos_diagnostics"] = Vector{Dict{Any, Any}}())
        push!(
            config_dict["extra_atmos_diagnostics"],
            Dict("short_name" => ["toa_fluxes_net"], "reduction_time" => "average", "period" => period),
        )
    end
    config_dict["calendar_dt"] = calendar_dt
    return nothing
end
