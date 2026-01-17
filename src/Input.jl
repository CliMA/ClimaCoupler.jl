"""
    Input

Module for parsing command-line arguments and configuration files.
"""
module Input

import ArgParse
import YAML
import Dates
import ClimaAtmos as CA
import ClimaUtilities.TimeManager: ITime
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaCore as CC
import ClimaCoupler
import ..Interfacer
import ..Utilities

export argparse_settings,
    parse_commandline, get_coupler_config_dict, get_coupler_args, get_land_fraction

const MODE_NAME_DICT = Dict(
    "amip" => Interfacer.AMIPMode,
    "cmip" => Interfacer.CMIPMode,
    "slabplanet" => Interfacer.SlabplanetMode,
    "slabplanet_aqua" => Interfacer.SlabplanetAquaMode,
    "slabplanet_terra" => Interfacer.SlabplanetTerraMode,
    "subseasonal" => Interfacer.SubseasonalMode,
)

"""
    argparse_settings()

Create and return an `ArgParseSettings` object with all command-line arguments
for ClimaCoupler simulations. Each option should include an argument type,
a default value, and a brief help string including the valid values for this option.
"""
function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        ### ClimaCoupler flags
        # Simulation-identifying information
        "--config_file"
        help = "A yaml file used to set the configuration of the coupled model [\"config/ci_configs/amip_default.yml\" (default)]"
        arg_type = String
        default = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml")
        "--job_id"
        help = "A unique identifier for this run, defaults to the config file name"
        arg_type = String
        default = nothing
        "--print_config_dict"
        help = "Boolean flag indicating whether to print the final configuration dictionary [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--mode_name"
        help = "Mode of coupled simulation. [`cmip`, `amip` (default), `subseasonal`, `slabplanet`, `slabplanet_aqua`, `slabplanet_terra`]"
        arg_type = String
        default = "amip"
        "--coupler_toml"
        help = "An optional list of paths to toml files used to overwrite the default model parameters."
        arg_type = Vector{String}
        default = []
        # Computational simulation setup information
        "--unique_seed"
        help = "Boolean flag indicating whether to set the random number seed to a unique value [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--FLOAT_TYPE"
        help = "Floating point precision  [`Float64` (default), `Float32`]"
        arg_type = String
        default = "Float64"
        "--device"
        help = "Device type to use [\"auto\" (default), \"CPUSingleThreaded\", \"CPUMultiThreaded\", \"CUDADevice\"]"
        arg_type = String
        default = "auto"
        # Time information
        "--use_itime"
        help = "Boolean flag indicating whether to use ITime (integer time) or not (will use Float64) [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--t_end"
        help = "End time of the simulation [\"800secs\"; allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        default = "800secs"
        "--t_start"
        help = "Start time of the simulation [\"0secs\" (default); allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        default = "0secs"
        "--start_date"
        help = "Start date of the simulation, in format \"YYYYMMDD\" [\"20100101\" (default)]"
        arg_type = String
        default = "20000101"
        "--dt_cpl"
        help = "Coupling time step in seconds [400 (default); allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        default = "400secs"
        "--dt"
        help = "Component model time step [allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        default = "400secs"
        "--dt_atmos"
        help = "Atmos simulation time step (alternative to `dt`; no default) [allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        "--dt_land"
        help = "Land simulation time step (alternative to `dt`; no default) [allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        "--dt_ocean"
        help = "Ocean simulation time step (alternative to `dt`; no default) [allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        "--dt_seaice"
        help = "Sea ice simulation time step (alternative to `dt`; no default) [allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        "--checkpoint_dt"
        help = "Time interval for checkpointing [\"90days\" (default)]"
        arg_type = String
        default = "90days"
        # Space information
        "--h_elem"
        help = "Number of horizontal elements to use for the boundary space [16 (default)]"
        arg_type = Int
        default = 16
        "--share_surface_space"
        help = "Boolean flag indicating whether to share the surface space between the surface models, atmosphere, and boundary [`true` (default), `false`]"
        arg_type = Bool
        default = true
        # Restart information
        "--detect_restart_files"
        help = "Boolean flag indicating whether to automatically use restart files if available [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--restart_dir"
        help = "Directory containing restart files"
        arg_type = String
        default = nothing
        "--restart_t"
        help = "Time in seconds rounded to the nearest index to use at `t_start` for restarted simulation [nothing (default)]"
        arg_type = Int
        default = nothing
        "--restart_cache"
        help = "Boolean flag indicating whether to read the cache from the restart file if available [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--save_cache"
        help = "Boolean flag indicating whether to save the state and cache or only the state when checkpointing [`true` (default), `false`]"
        arg_type = Bool
        default = true
        # Diagnostics information
        "--use_coupler_diagnostics"
        help = "Boolean flag indicating whether to compute and output coupler diagnostics [`true` (default), `false`]"
        arg_type = Bool
        default = true
        # Physical simulation information
        "--evolving_ocean"
        help = "Boolean flag indicating whether to use a dynamic slab ocean model, as opposed to constant surface temperatures [`true` (default), `false`]"
        arg_type = Bool
        default = true
        # Conservation and RMSE check information
        "--energy_check"
        help = "Boolean flag indicating whether to check energy conservation [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--conservation_softfail"
        help = "Boolean flag indicating whether to soft fail on conservation errors [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--rmse_check"
        help = "Boolean flag indicating whether to check RMSE of some physical fields [`false` (default), `true`]"
        arg_type = Bool
        default = false
        # Output information
        "--coupler_output_dir"
        help = "Directory to save output files. Note that TempestRemap fails if interactive and paths are too long. [\"output\" (default)]"
        arg_type = String
        default = "output"
        # ClimaAtmos specific
        "--surface_setup"
        help = "Triggers ClimaAtmos into the coupled mode [`PrescribedSurface` (default), `DefaultMoninObukhov`]" # retained here for standalone Atmos benchmarks
        arg_type = String
        default = "PrescribedSurface"
        "--atmos_config_file"
        help = "An optional YAML file used to overwrite the default model parameters."
        arg_type = String
        default = nothing
        "--atmos_log_progress"
        help = "Use the ClimaAtmos walltime logging callback instead of the default ClimaCoupler one [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--albedo_model"
        help = "Type of albedo model. [`ConstantAlbedo`, `RegressionFunctionAlbedo`, `CouplerAlbedo` (default)]"
        arg_type = String
        default = "CouplerAlbedo"
        "--extra_atmos_diagnostics"
        help = "List of dictionaries containing information about additional atmosphere diagnostics to output [nothing (default)]"
        arg_type = Vector{Dict{Any, Any}}
        default = []
        ### ClimaLand specific
        "--land_model"
        help = "Land model to use. [`bucket` (default), `integrated`]"
        arg_type = String
        default = "bucket"
        "--land_temperature_anomaly"
        help = "Type of temperature anomaly for land model. [`amip`, `aquaplanet` (default), `nothing`]"
        arg_type = String
        default = "aquaplanet"
        "--use_land_diagnostics"
        help = "Boolean flag indicating whether to compute and output land model diagnostics [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--land_spun_up_ic"
        help = "Boolean flag to indicate whether to use integrated land initial conditions from spun up state [`true` (default), `false`]"
        arg_type = Bool
        default = true
        # BucketModel specific
        "--bucket_albedo_type"
        help = "Access bucket surface albedo information from data file. [`map_static` (default), `function`, `map_temporal`, `era5`]"
        arg_type = String
        default = "map_static" # to be replaced by land config file, when available
        "--bucket_initial_condition"
        help = "A file path for a NetCDF file (read documentation about requirements)"
        arg_type = String
        default = ""
        "--era5_initial_condition_dir"
        help = "Directory containing ERA5 initial condition files (subseasonal mode). Filenames inferred from start_date [none (default)]. Generated with `https://github.com/CliMA/WeatherQuest`"
        arg_type = String
        default = nothing
        # Ice model specific
        "--ice_model"
        help = "Sea ice model to use. [`prescribed` (default), `clima_seaice`]"
        arg_type = String
        default = "prescribed"
        "--land_fraction_source"
        help = "Source for land fraction data. [`etopo` (default) uses ETOPO-derived landsea_mask artifact, `era5` uses ERA5 land fraction artifact]"
        arg_type = String
        default = "etopo"
        "--binary_area_fraction"
        help = "Boolean flag indicating whether to use binary (thresholded) area fractions for land and ice [`true` (default), `false`]. When true, land fraction > eps becomes 1, and ice fraction > 0.5 becomes 1."
        arg_type = Bool
        default = true
    end
    return s
end

"""
    parse_commandline(settings)

Parse command-line arguments using the provided `ArgParseSettings` object.

# Arguments
- `settings`: An `ArgParse.ArgParseSettings` object containing the command-line arguments to parse

# Returns
- A dictionary of parsed command-line arguments
"""
parse_commandline(settings) = ArgParse.parse_args(ARGS, settings)

"""
    get_coupler_config_dict(config_file)

Read in the configuration file and job ID from the command line.
A dictionary is constructed from the input configuration file and returned.

Since the atmosphere model also uses a configuration file, we read in the atmosphere
configuration file specified in the coupler configuration file (if any), and overwrite
it with the coupler configuration.

The order of priority for overwriting configuration options from lowest to highest is:
    1. ClimaAtmos defaults
    2. ClimaCoupler defaults (defined in `Input.argparse_settings()`)
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
        isnothing(job_id) ? splitext(basename(config_file))[1] : job_id

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

# Arguments
- `config_dict`: A dictionary mapping configuration keys to the specified settings

# Returns
- A NamedTuple of all arguments needed for the coupled simulation
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
    sim_mode = MODE_NAME_DICT[mode_name]
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

    # Land fraction source
    land_fraction_source = config_dict["land_fraction_source"]

    # Binary area fraction
    binary_area_fraction = config_dict["binary_area_fraction"]

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
        land_fraction_source,
        binary_area_fraction,
    )
end

### Helper functions used in argument parsing ###

"""
    get_diag_period(t_start, t_end)

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
    if all(
        key -> haskey(config_dict, key) && !isnothing(config_dict[key]),
        component_dt_names,
    )
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
            if haskey(config_dict, key) && !isnothing(config_dict[key])
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
    get_land_fraction(boundary_space, comms_ctx; land_fraction_source = "etopo", binary_area_fraction = true)

Read and remap the land-sea fraction field onto the coupler boundary grid.

# Arguments
- `boundary_space`: The boundary space onto which to remap the land fraction.
- `comms_ctx`: The communications context.
- `land_fraction_source`: Source of land fraction data. Either "etopo" (default) or "era5".
- `binary_area_fraction`: If true (default), threshold land fraction to binary (0 or 1).

# Returns
- A field containing land fraction values (0 to 1) on the boundary space.

Note: 
Land-sea Fraction
    This is a static field that contains the area fraction of land and sea, ranging from 0 to 1.
    If applicable, sea ice is included in the sea fraction at this stage.
    Note that land-sea area fraction is different to the land-sea mask, which is a binary field
    (masks are used internally by the coupler to indicate passive cells that are not populated by a given component model).

    Two sources are supported via the `land_fraction_source` config option:
    - "etopo": ETOPO-derived binary land-sea mask (landsea_mask_60arcseconds artifact)
    - "era5": ERA5 land fraction field (era5_land_fraction artifact)
"""
function get_land_fraction(
    boundary_space,
    comms_ctx;
    land_fraction_source::String = "etopo",
    binary_area_fraction::Bool = true,
)
    FT = CC.Spaces.undertype(boundary_space)

    if land_fraction_source == "era5"
        land_fraction_data = joinpath(
            @clima_artifact("era5_land_fraction", comms_ctx),
            "era5_land_fraction.nc",
        )
        land_fraction = SpaceVaryingInput(land_fraction_data, "lsm", boundary_space)
    elseif land_fraction_source == "etopo"
        land_fraction_data = joinpath(
            @clima_artifact("landsea_mask_60arcseconds", comms_ctx),
            "landsea_mask.nc",
        )
        land_fraction = SpaceVaryingInput(land_fraction_data, "landsea", boundary_space)
    else
        error(
            "Unknown land_fraction_source: $land_fraction_source. Must be \"etopo\" or \"era5\".",
        )
    end

    # Ensure land fraction is finite/not NaN and clamp to [0, 1]
    land_fraction = ifelse.(isfinite.(land_fraction), land_fraction, FT(0))
    land_fraction = max.(min.(land_fraction, FT(1)), FT(0))

    if binary_area_fraction
        land_fraction = ifelse.(land_fraction .> eps(FT), FT(1), FT(0))
    else
        land_fraction = ifelse.(land_fraction .> eps(FT), land_fraction, FT(0))
    end

    return land_fraction
end

end # module Input
