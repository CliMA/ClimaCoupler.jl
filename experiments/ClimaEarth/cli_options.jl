import ArgParse
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
        help = "Directory to save output files. Note that TempestRemap fails if interactive and paths are too long. [\"experiments/ClimaEarth/output\" (default)]"
        arg_type = String
        default = "experiments/ClimaEarth/output"
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
        help = "Access bucket surface albedo information from data file. [`map_static` (default), `function`, `map_temporal`]"
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
    end
    return s
end

parse_commandline(s) = ArgParse.parse_args(ARGS, s)
