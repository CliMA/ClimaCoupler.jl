import ArgParse
function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        ### ClimaCoupler flags
        # Simulation-identifying information
        "--config_file"
        help = "A yaml file used to set the configuration of the coupled model [\"config/ci_configs/amip_default.yml\" (default)]"
        arg_type = String
        default = "config/ci_configs/amip_default.yml"
        "--job_id"
        help = "A unique identifier for this run, defaults to the config file name"
        arg_type = String
        default = nothing
        "--print_config_dict"
        help = "Boolean flag indicating whether to print the final configuration dictionary [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--mode_name"
        help = "Mode of coupled simulation. [`amip` (default), `slabplanet`, `slabplanet_aqua`, `slabplanet_terra`, `slabplanet_eisenman`]"
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
        "--dt_save_to_sol"
        help = "Time interval for saving output [\"10days\" (default); allowed formats: \"Nsecs\", \"Nmins\", \"Nhours\", \"Ndays\", \"Inf\"]"
        arg_type = String
        default = "10days"
        "--checkpoint_dt"
        help = "Time interval for hourly checkpointing [\"20days\" (default)]"
        arg_type = String
        default = "20days"
        # Restart information
        "--restart_dir"
        help = "Directory containing restart files"
        arg_type = String
        default = nothing
        "--restart_t"
        help = "Time in seconds rounded to the nearest index to use at `t_start` for restarted simulation [0 (default)]"
        arg_type = Int
        default = 0
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
        "--mono_surface"
        help = "Boolean flag indicating whether (1st order) monotone and conservative remapping is applied. [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--turb_flux_partition"
        help = "Method to partition turbulent fluxes. [`PartitionedStateFluxes`, `CombinedStateFluxesMOST` (default)]"
        arg_type = String
        default = "CombinedStateFluxesMOST"
        # Conservation information
        "--energy_check"
        help = "Boolean flag indicating whether to check energy conservation [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--conservation_softfail"
        help = "Boolean flag indicating whether to soft fail on conservation errors [`false` (default), `true`]"
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
        "--albedo_model"
        help = "Type of albedo model. [`ConstantAlbedo`, `RegressionFunctionAlbedo`, `CouplerAlbedo` (default)]"
        arg_type = String
        default = "CouplerAlbedo"
        "--atmos_config_repo"
        help = "The repository containing the ClimaAtmos configuration file to use [`ClimaAtmos` (default), `ClimaCoupler`]"
        arg_type = String
        default = "ClimaAtmos"
        "--extra_atmos_diagnostics"
        help = "List of dictionaries containing information about additional atmosphere diagnostics to output [nothing (default)]"
        arg_type = Vector{Dict{Any, Any}}
        default = []
        ### ClimaLand specific
        "--land_domain_type"
        help = "Type of land domain. [`sphere` (default), `single_column`]"
        arg_type = String
        default = "sphere"
        "--land_albedo_type"
        help = "Access land surface albedo information from data file. [`map_static` (default), `function`, `map_temporal`]"
        arg_type = String
        default = "map_static" # to be replaced by land config file, when available
        "--land_initial_condition"
        help = "A file path for a NetCDF file (read documentation about requirements)"
        arg_type = String
        default = ""
        "--land_temperature_anomaly"
        help = "Type of temperature anomaly for bucket model. [`amip`, `aquaplanet` (default)]"
        arg_type = String
        default = "aquaplanet"
        "--use_land_diagnostics"
        help = "Boolean flag indicating whether to compute and output land model diagnostics [`true` (default), `false`]"
        arg_type = Bool
        default = true
    end
    return s
end

parse_commandline(s) = ArgParse.parse_args(ARGS, s)
