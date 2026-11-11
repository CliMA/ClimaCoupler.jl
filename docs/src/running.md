# Running a Simulation

This page walks through the most common ways to set up and run a simulation with
ClimaCoupler.jl. For a full description of all available configuration options, see the
[Input](@ref) module documentation.

The setup for other simulation types (slabplanet, CMIP, etc) is the same.
The only pieces that need to change are the selected configuration file and the
Julia environment. Each environment has a `code_loading.jl` script that loads
ClimaCoupler together with the packages that trigger its extensions:

| Simulation type     | Julia environment   | Package loading script             |
| ------------------- | ------------------- | ---------------------------------- |
| AMIP, slabplanet    | `experiments/AMIP/` | `experiments/AMIP/code_loading.jl` |
| CMIP                | `experiments/CMIP/` | `experiments/CMIP/code_loading.jl` |

The CMIP script additionally loads Oceananigans, ClimaSeaIce, KernelAbstractions,
ConservativeRegridding, and Adapt, which trigger the CMIP extensions.

## Setting up the environment

Clone the repository and run all commands on this page from its root directory: config
file paths and the `output/` directory are resolved relative to the directory Julia is
started from.

```bash
git clone https://github.com/CliMA/ClimaCoupler.jl.git
cd ClimaCoupler.jl
```

Use Julia 1.10 or 1.13. Each experiment environment ships manifests pinned to those versions
(`Manifest.toml` for 1.10, `Manifest-v1.13.toml` for 1.13); other Julia versions fall back
to the 1.10 manifest and may need to re-resolve it. With
[juliaup](https://github.com/JuliaLang/juliaup), `juliaup add 1.13` installs a matching
version and `julia +1.13` selects it.

Install the packages of the environment you will use, once and again after pulling changes
to its manifest:

```bash
julia --project=experiments/AMIP -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/CMIP -e 'using Pkg; Pkg.instantiate()'  # only for CMIP
```

The first simulation also downloads its input data (land-sea masks, SST and sea-ice data,
etc.) as Julia artifacts.

## From the command line

The simplest way to launch an AMIP simulation is from a shell in the root of the
ClimaCoupler.jl repository:

```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl
```

This uses the default configuration at `config/ci_configs/amip_default.yml`. To use a
different config file, pass the relative path via `--config_file`:

```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl --config_file="config/ci_configs/amip_default.yml"
```

CMIP simulations use the CMIP environment and driver, and must be given a CMIP config
file explicitly, since the default config is an AMIP setup:

```bash
julia --project=experiments/CMIP experiments/CMIP/run_simulation.jl --config_file="config/ci_configs/cmip_oceananigans_climaseaice.yml"
```

A collection of ready-to-use config files for common setups can be found in
`config/ci_configs/`. Note that many of these simulations have a very short
run time as they're primarily used for software testing. This can be changed by editing
the files directly, or by modifying the produced config dictionary as described below.

A list of all available configuration options and their defaults can be found on the
[Input](@ref) page.

## From the REPL

To run a simulation interactively, start Julia using the AMIP project
and include the package loading script, which imports all required packages and triggers
the extensions needed for AMIP and slabplanet runs:

!!! note "Working directory"
    All REPL examples below assume Julia was started from the **repository root**
    directory, e.g. `julia --project=experiments/AMIP`.

!!! note "CMIP"
    For CMIP simulations, start Julia with `julia --project=experiments/CMIP`,
    replace `experiments/AMIP/code_loading.jl` with `experiments/CMIP/code_loading.jl`,
    and pass a CMIP config file such as
    `config/ci_configs/cmip_oceananigans_climaseaice.yml` to `CoupledSimulation`.

```julia
include("experiments/AMIP/code_loading.jl")

cs = CoupledSimulation()  # uses amip_default.yml
run!(cs)
postprocess(cs)
```

To use a specific config file:

```julia
include("experiments/AMIP/code_loading.jl")

config_file = "config/ci_configs/amip_default.yml"
cs = CoupledSimulation(config_file)
run!(cs)
postprocess(cs)
```

## Running on GPU or with MPI

Both experiment environments include CUDA.jl and MPI.jl, and `code_loading.jl` loads them
through `ClimaComms.@import_required_backends` when these environment variables request them:

```bash
export CLIMACOMMS_DEVICE="CUDA"   # run on an NVIDIA GPU
export CLIMACOMMS_CONTEXT="MPI"   # run distributed, with one Julia process per MPI rank
```

The macro reads the variables when it runs, so in the REPL set them with `ENV` *before*
including `code_loading.jl`:

```julia
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
include("experiments/AMIP/code_loading.jl")
```

Keep the `device` configuration option at its default, `"auto"`, so that
`CLIMACOMMS_DEVICE` selects the device. If a single-process run crashes with an MPI error
(common in interactive SLURM sessions), set `CLIMACOMMS_CONTEXT="SINGLETON"`. For CUDA
driver compatibility and how to verify that a run is on the GPU, see
[Running a CliMA Model on GPU](https://github.com/CliMA/DeveloperGuides/blob/main/workflow/running_on_gpu.md).

## Modifying configuration interactively

For cases where you want to adjust a few settings without creating a new YAML file, you
can load a config into a dictionary, modify it, and pass it directly to
`CoupledSimulation`. For example, to run for one day instead of the default:

```julia
include("experiments/AMIP/code_loading.jl")

config_file = "config/ci_configs/amip_default.yml"
config_dict = Input.get_coupler_config_dict(config_file)
config_dict["t_end"] = "1days"

cs = CoupledSimulation(config_dict)
run!(cs)
postprocess(cs)
```

Any key from the [Available Configuration Options](@ref) table can be modified this way.
The dictionary uses the same key names as the YAML config file and the CLI flags.

## Using ClimaCoupler as an installed package

If you are using ClimaCoupler.jl as a package dependency rather than working inside a
clone of the repository, you just need to load ClimaCoupler and trigger the required
extensions — the equivalent of what `code_loading.jl` does. These trigger packages are
weak dependencies of ClimaCoupler, so they must be added to your own environment
(e.g. with `Pkg.add`):

```julia
using ClimaCoupler

# Trigger ClimaCouplerMakieExt (plotting)
using CairoMakie, GeoMakie, Makie, Poppler_jll, Printf

# Trigger ClimaCouplerCMIPExt (only needed for CMIP simulations).
# Together with the Makie packages above, Oceananigans also triggers ClimaCouplerCMIPMakieExt.
import Oceananigans, ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt

# Trigger ClimaCouplerClimaLandExt
import ClimaLand

# Trigger ClimaCouplerClimaAtmosExt (also required by ClimaCouplerCMIPExt)
import ClimaAtmos

# Load the backends (e.g. CUDA, MPI) required by the current ClimaComms context
import ClimaComms
ClimaComms.@import_required_backends
```

Once the packages are loaded, construct and run the simulation using a config file:

```julia
config_file = joinpath(pkgdir(ClimaCoupler), "config", "ci_configs", "amip_default.yml")

cs = CoupledSimulation(config_file)
run!(cs)
postprocess(cs)
```

## Stepping manually

Instead of `run!`, you can drive the simulation step by step using `step!`. Each call
advances the simulation by one coupling timestep (`cs.Δt_cpl`), which lets you inspect
or modify state between steps:

```julia
include("experiments/AMIP/code_loading.jl")

cs = CoupledSimulation()

# Advance the simulation by three timesteps
step!(cs)
@info "current simulation time: $(cs.t[])"
step!(cs)
@info "current simulation time: $(cs.t[])"
step!(cs)
@info "current simulation time: $(cs.t[])"
```

To run the simulation for its entire duration using individual steps, you can
do the following:

```julia
include("experiments/AMIP/code_loading.jl")

cs = CoupledSimulation()

while cs.t[] < cs.tspan[end]
    step!(cs)
end
```

This is equivalent to what `run!` does internally, without the precompilation warmup
and SYPD timing. See [SimCoordinator](@ref) for the full API.

## Single-column model (SCM) experiments

To run a single-column model (SCM) experiment, set `domain_type: "column"`, **`scm_surface_type`**
(which surface runs: land, ocean, or sea ice), and **`column_latlon`** (where the column sits).
Example config file:

```yaml
domain_type: "column"
column_latlon: [37.0, -120.0]   # [latitude, longitude] in degrees
scm_surface_type: "ocean"       # "land", "ocean", or "sea_ice"
mode_name: "slabplanet_aqua"    # or "amip", "slabplanet_terra", etc.
dt: "120secs"
dt_cpl: "120secs"
t_end: "480secs"
start_date: "20100101"
# ... other options (rad, microphysics_model, surface_setup, vert_diff, z_elem, etc.)
```

Run it the same way as a global run, e.g.:

```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl --config_file=config/ci_configs/scm_slabplanet_aqua.yml --job_id=my_scm
```

Or from the REPL:

```julia
include("experiments/AMIP/code_loading.jl")
cs = CoupledSimulation("config/ci_configs/scm_slabplanet_aqua.yml")
run!(cs)
```

Ready-to-use SCM configs are in `config/ci_configs/` (e.g. `scm_slabplanet_aqua.yml`,
`scm_slabplanet_terra_land.yml`, `scm_amip_ocean.yml`).

### SCM surface model selection

`scm_surface_type` is used to select which surface model runs in SCM mode. It must
be set to `"land"`, `"ocean"`, or `"sea_ice"`. The coupler then keeps exactly one of the land,
ocean, or sea-ice component models active (the others are disabled). Which *kind* of ocean,
land, or ice model (e.g. slab vs prescribed) still follows from `mode_name` and the usual
`ocean_model` / `land_model` / `ice_model` config where applicable. Multiple surface
model types (e.g. half ocean, half land) are currently not supported in SCM mode.

`column_latlon` sets the geographic location of the column (e.g. for spatially-varying
land parameters, prescribed ocean albedo calculation, radiation).
