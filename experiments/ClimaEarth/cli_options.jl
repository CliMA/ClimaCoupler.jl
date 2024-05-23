import ArgParse
function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        # ClimaCoupler flags
        "--dt_cpl"
        help = " Coupling time step in seconds"
        arg_type = Int
        default = 400
        "--anim"
        help = "Boolean flag indicating whether to make animations"
        arg_type = Bool
        default = false
        "--energy_check"
        help = "Boolean flag indicating whether to check energy conservation"
        arg_type = Bool
        default = false
        "--ci_plots"
        help = "Boolean flag indicating whether to make CI plots"
        arg_type = Bool
        default = false
        "--conservation_softfail"
        help = "Boolean flag indicating whether to soft fail on conservation errors"
        arg_type = Bool
        default = false
        "--mode_name"
        help = "Mode of coupled simulation. [`amip`, `slabplanet`, `slabplanet_aqua`, `slabplanet_terra`, `slabplanet_eisenman`]"
        arg_type = String
        default = "amip"
        "--mono_surface"
        help = "Boolean flag indicating whether (1st order) monotone and conservative remapping is applied."
        arg_type = Bool
        default = false
        "--turb_flux_partition"
        help = "Method to partition turbulent fluxes. [`PartitionedStateFluxes`, `CombinedStateFluxesMOST`]"
        arg_type = String
        default = "CombinedStateFluxesMOST"
        "--hourly_checkpoint"
        help = "Boolean flag indicating whether to checkpoint at intervals of 1 hour or multiple hours"
        arg_type = Bool
        default = false
        "--hourly_checkpoint_dt"
        help = "Time interval for hourly checkpointing in hours (20 days by default)"
        arg_type = Int
        default = 480
        "--coupler_output_dir"
        help = "Directory to save output files. Note that TempestRemap fails if interactive and paths are too long."
        arg_type = String
        default = "experiments/ClimaEarth/output"
        "--restart_dir"
        help = "Directory containing restart files"
        arg_type = String
        default = "unspecified"
        "--restart_t"
        help = "Restart time"
        arg_type = Int
        default = 0
        "--config_file"
        help = "A yaml file used to set the configuration of the coupled model"
        arg_type = String
        "--job_id"
        help = "A unique identifier for this run"
        arg_type = String
        default = nothing
        "--print_config_dict"
        help = "Boolean flag indicating whether to print the final configuration dictionary"
        arg_type = Bool
        default = true
        "--FLOAT_TYPE"
        help = "Floating point precision  [`Float64` (default), `Float32`]"
        arg_type = String
        default = "Float64"
        "--coupler_toml_file"
        help = "A toml file used to overwrite the model parameters. If nothing is specified, the default parameters are used."
        "--evolving_ocean"
        help = "Boolean flag indicating whether to use a dynamic slab ocean model or constant surface temperatures"
        arg_type = Bool
        default = true
        "--device"
        help = "Device type to use [`auto` (default) `CPUSingleThreaded`, `CPUMultiThreaded`, `CUDADevice`]"
        arg_type = String
        default = "auto"
        "--use_coupler_diagnostics"
        help = "Boolean flag indicating whether to compute and output coupler diagnostics [`true` (default), `false`]"
        arg_type = Bool
        default = true
        # ClimaAtmos specific
        "--surface_setup"
        help = "Triggers ClimaAtmos into the coupled mode [`PrescribedSurface` (default)]" # retained here for standalone Atmos benchmarks
        arg_type = String
        default = "PrescribedSurface"
        "--atmos_config_file"
        help = "A yaml file used to set the atmospheric model configuration. If nothing is specified, the default configuration is used."
        "--albedo_model"
        help = "Type of albedo model. [`ConstantAlbedo` (default), `RegressionFunctionAlbedo`, `CouplerAlbedo`]"
        arg_type = String
        default = "CouplerAlbedo"
        "--atmos_config_repo"
        help = "The repository containing the ClimaAtmos configuration file to use [`ClimaAtmos` (default), `ClimaCoupler`]"
        arg_type = String
        default = "ClimaAtmos"
        # ClimaLand specific
        "--land_albedo_type"
        help = "Access land surface albedo information from data file. [`function`, `map_static`, `map_temporal`]"
        arg_type = String
        default = "map_static" # to be replaced by land config file, when available
        "--land_domain_type"
        help = "Type of land domain. [`sphere` (default), `single_column`]"
        arg_type = String
        default = "sphere"
        "--land_temperature_anomaly"
        help = "Type of temperature anomaly for bucket model. [`amip`, `aquaplanet` (default)]"
        arg_type = String
        default = "aquaplanet"
    end
    return s
end

parse_commandline(s) = ArgParse.parse_args(ARGS, s)
