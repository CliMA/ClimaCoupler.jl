# Input

The `Input` module provides functions for parsing command-line arguments and loading
configuration files for ClimaCoupler simulations.

## Providing Configuration Options

ClimaCoupler simulations have a variety of configuration options, which are
described in the section [Available Configuration Options](#available-configuration-options) below.

Users can configure a coupled simulation through 2 main entry points:
command-line arguments and configuration files. Each method, and the
priority between them, is explained here.

Note that all of these configuration are optional; a default AMIP simulation
can be set up with the 0-argument constructor `CoupledSimulation()`.

### Configuration Files

Configuration files are YAML files that specify simulation parameters. They can be
used to specify a valid setting for any configuration option.

Note that since the atmosphere model has many configuration options on its own,
we also support providing an atmosphere-specific configuration file. As a result,
multiple configuration files can be used together:

- **Coupler config file** (`--config_file`): Main configuration file for the coupled simulation
- **Atmos config file** (`--atmos_config_file`): Optional ClimaAtmos-specific configuration
- **TOML parameter files** (`--coupler_toml`): One or more TOML files containing model parameters

When multiple config files are specified, values in the coupler config file will take
precedence over those in the atmosphere config file. This is explained in more detail
in the [Precendence of Config Files and CLI Arguments](#precendence-of-config-files-and-cli-arguments)
section below.

### Parameter TOML Files

Analogous to configuration files, which specify simulation options such as component models or
parameterizations to use, ClimaCoupler also accepts optional TOML files of parameter values.

If a TOML file is provided, the default parameter values will be overwritten by the values
specified in the TOML. If no TOML file is provided, default parameter values will be used.

Since ClimaCoupler accepts coupler config files _and_ atmosphere config files, we may encounter
a case where a TOML file is specified in both config files. In this case, the coupler TOML
file takes highest priority; only if there is no coupler TOML will the atmosphere-specific
TOML be used.

### Command-Line Input (CLI) Arguments

Command-line arguments can be provided when running a simulation:

```bash
julia run_amip.jl --config_file="path/to/config.yml" --job_id="amip_default"
```

Typically, we rely mostly on configuration files, and provide only the `config_file`
and `job_id` via CLI arguments, though if desired any input can be specified in the command line.
All available options and their defaults can be viewed by running with `--help`.

### Precendence of Config Files and CLI Arguments

Users have the ability to specify arguments via both CLI arguments and configuration
files for any given simulation. If this is done, the options are merged
with the following precedence (from lowest to highest priority):

1. **ClimaAtmos defaults** - Default values from the ClimaAtmos package
2. **ClimaCoupler defaults** - Default values defined in [`argparse_settings()`](@ref)
3. **Command-line arguments** - Arguments passed via the command line
4. **ClimaAtmos configuration file** - YAML file specified via `--atmos_config_file` or in the coupler config file
5. **ClimaCoupler configuration file** - YAML file specified via `--config_file` (default: `config/ci_configs/amip_default.yml`)

## Available Configuration Options

The following table lists all available command-line arguments organized by category:

#### Simulation-identifying information

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--config_file` | String | `config/ci_configs/amip_default.yml` | Any valid file path | YAML file used to set the configuration of the coupled model |
| `--job_id` | String | `nothing` | Any string | A unique identifier for this run, defaults to the config file name |
| `--print_config_dict` | Bool | `true` | `true`, `false` | Whether to print the final configuration dictionary |
| `--mode_name` | String | `"amip"` | `cmip`, `amip`, `subseasonal`, `slabplanet`, `slabplanet_aqua`, `slabplanet_terra` | Mode of coupled simulation |
| `--coupler_toml` | Vector{String} | `[]` | Any list of valid TOML file paths | Optional list of paths to TOML files used to overwrite default model parameters |

#### Computational simulation setup

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--unique_seed` | Bool | `false` | `true`, `false` | Whether to set the random number seed to a unique value |
| `--FLOAT_TYPE` | String | `"Float64"` | `Float64`, `Float32` | Floating point precision |
| `--device` | String | `"auto"` | `auto`, `CPUSingleThreaded`, `CPUMultiThreaded`, `CUDADevice` | Device type to control running on CPU or GPU |

#### Time information

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--use_itime` | Bool | `true` | `true`, `false` | Whether to use ClimaUtilities ITime (integer time) or Float64 |
| `--t_end` | String | `"800secs"` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | End time of the simulation, relative to the start date |
| `--t_start` | String | `"0secs"` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | Start time of the simulation, relative to the start date |
| `--start_date` | String | `"20000101"` | `"YYYYMMDD"` format | Start date of the simulation |
| `--dt_cpl` | String | `"400secs"` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | Coupling time step |
| `--dt` | String | `"400secs"` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | Component model time step (used if individual component dt's not specified) |
| `--dt_atmos` | String | `nothing` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | Atmos simulation time step (alternative to `dt`) |
| `--dt_land` | String | `nothing` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | Land simulation time step (alternative to `dt`) |
| `--dt_ocean` | String | `nothing` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | Ocean simulation time step (alternative to `dt`) |
| `--dt_seaice` | String | `nothing` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | Sea ice simulation time step (alternative to `dt`) |
| `--checkpoint_dt` | String | `"90days"` | `"Nsecs"`, `"Nmins"`, `"Nhours"`, `"Ndays"`, `"Inf"` | Time interval for checkpointing |

Note: If any component model-specific timestep is specified, _all_ component-model
specific timesteps should be specified, rather than only `dt`.

#### Space information

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--h_elem` | Int | `16` | Any positive integer | Number of horizontal elements to use for the boundary space |
| `--share_surface_space` | Bool | `true` | `true`, `false` | Whether to share the surface space between surface models, atmosphere, and boundary |

#### Restart information

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--detect_restart_files` | Bool | `false` | `true`, `false` | Whether to automatically use restart files if available |
| `--restart_dir` | String | `nothing` | Any valid directory path | Directory containing restart files |
| `--restart_t` | Int | `nothing` | Any integer (seconds) | Time in seconds rounded to nearest index to use at `t_start` for restarted simulation |
| `--restart_cache` | Bool | `true` | `true`, `false` | Whether to read the cache from the restart file if available |
| `--save_cache` | Bool | `true` | `true`, `false` | Whether to save the state and cache or only the state when checkpointing |

#### Diagnostics and output

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--use_coupler_diagnostics` | Bool | `true` | `true`, `false` | Whether to compute and output coupler diagnostics |
| `--coupler_output_dir` | String | `"experiments/ClimaEarth/output"` | Any valid directory path | Directory to save output files |



#### Conservation and RMSE checks

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--energy_check` | Bool | `false` | `true`, `false` | Whether to check energy conservation |
| `--conservation_softfail` | Bool | `false` | `true`, `false` | Whether to soft fail on conservation errors |
| `--rmse_check` | Bool | `false` | `true`, `false` | Whether to check RMSE of some physical fields |

#### ClimaAtmos specific

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--surface_setup` | String | `"PrescribedSurface"` | `PrescribedSurface`, `DefaultMoninObukhov` | Triggers ClimaAtmos into coupled mode |
| `--atmos_config_file` | String | `nothing` | Any valid file path | Optional YAML file used to overwrite default model parameters |
| `--atmos_log_progress` | Bool | `false` | `true`, `false` | Use ClimaAtmos walltime logging callback instead of default ClimaCoupler one |
| `--albedo_model` | String | `"CouplerAlbedo"` | `ConstantAlbedo`, `RegressionFunctionAlbedo`, `CouplerAlbedo` | Type of albedo model |
| `--extra_atmos_diagnostics` | Vector{Dict{Any, Any}} | `[]` | List of dictionaries | List of dictionaries containing information about additional atmosphere diagnostics to output |

#### ClimaLand specific

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--land_model` | String | `"bucket"` | `bucket`, `integrated` | Land model to use |
| `--land_temperature_anomaly` | String | `"aquaplanet"` | `amip`, `aquaplanet`, `nothing` | Type of temperature anomaly for land model |
| `--use_land_diagnostics` | Bool | `true` | `true`, `false` | Whether to compute and output land model diagnostics |
| `--land_spun_up_ic` | Bool | `true` | `true`, `false` | Whether to use integrated land initial conditions from spun up state |
| `--bucket_albedo_type` | String | `"map_static"` | `map_static`, `function`, `map_temporal`, `era5` | Access bucket surface albedo information from data file. Use `era5` for ERA5-derived processed albedo files (requires `era5_initial_condition_dir`) |
| `--bucket_initial_condition` | String | `""` | Any valid file path | File path for a NetCDF file (read documentation about requirements). In subseasonal mode, automatically inferred from `era5_initial_condition_dir` if not specified |
| `--era5_initial_condition_dir` | String | `nothing` | Any valid directory path | Directory containing ERA5 initial condition files (subseasonal mode). Filenames inferred from `start_date`. Generated with `https://github.com/CliMA/WeatherQuest` |
| `--land_fraction_source` | String | `"etopo"` | `etopo`, `era5` | Source for land fraction data. `etopo` uses ETOPO-derived landsea_mask artifact (binary), `era5` uses ERA5/IFS land fraction artifact (0.0 - 1.0), which includes large inland seas and lakes. |
| `--binary_area_fraction` | Bool | `true` | `true`, `false` | Whether to use binary (thresholded) area fractions for land and ice. When true, land fraction > eps becomes 1, and ice fraction > 0.5 becomes 1 |

#### Ice model specific

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--ice_model` | String | `"prescribed"` | `prescribed`, `clima_seaice` | Sea ice model to use |


#### Ocean model specific

| Argument | Type | Default | Valid Options | Description |
|----------|------|---------|---------------|-------------|
| `--evolving_ocean` | Bool | `true` | `true`, `false` | Whether to use a dynamic slab ocean model, as opposed to constant surface temperatures |

## Input API

```@docs
Input.argparse_settings
Input.parse_commandline
Input.get_coupler_config_dict
Input.get_coupler_args
Input.get_land_fraction
```
