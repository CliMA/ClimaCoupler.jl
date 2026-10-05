# # AMIP Driver

#=
## Overview

AMIP is a standard experimental protocol of the Program for Climate Model Diagnosis & Intercomparison (PCMDI).
It is used as a model benchmark for the atmospheric and land model components, while sea-surface temperatures (SST) and sea-ice concentration (SIC)
are prescribed using time-interpolations between monthly observed data. We use standard data files with original sources:
- SST and SIC: https://gdex.ucar.edu/dataset/158_asphilli.html
- land-sea mask: https://www.ncl.ucar.edu/Applications/Data/#cdf

For more information, see the PCMDI's specifications for [AMIP I](https://pcmdi.github.io/mips/amip/) and [AMIP II](https://pcmdi.github.io/mips/amip2/).

## Running the AMIP configuration
To run a coupled simulation in the default AMIP configuration, run the
following command from the root directory of the repository:
```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl
```

## Configuration
You can also specify a custom configuration file to run the coupled simulation
in a different setup. The configuration file should be a TOML file that overwrites
the input fields specified in the ClimaCoupler Input module.
A set of example configuration files can be found in the `config/ci_configs/` directory.

For example, to run the coupled simulation with a different configuration file:
```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl --config_file="path/to/config.toml"
```

To run the coupled simulation interactively with a different configuration file,
set the `config_file` variable in this script to be the path to that file.

For more details about running a coupled simulation, including how to run in a
Slabplanet configuration, please see our documentation.
=#

# Load the necessary modules to run the coupled simulation
include("code_loading.jl")

# The Makie plotting stack is heavy to load and compile, so by default it is
# deferred until after `run!` to keep the time-to-first-timestep low (see
# docs/src/precompilation_strategies.md). Opt out via environment variables:
#   CLIMACOUPLER_PLOTS_DURING_RUN=1  load Makie before the run (e.g. for plotting
#                                    callbacks that fire during the coupling loop)
#   CLIMACOUPLER_SKIP_POSTPROCESS=1  skip Makie and postprocessing entirely
plots_during_run = get(ENV, "CLIMACOUPLER_PLOTS_DURING_RUN", "0") in ("1", "true")
skip_postprocess = get(ENV, "CLIMACOUPLER_SKIP_POSTPROCESS", "0") in ("1", "true")
plots_during_run &&
    !skip_postprocess &&
    include(joinpath(@__DIR__, "..", "load_plotting.jl"))

# Get the configuration file from the command line (or manually set it here)
config_file = Input.parse_commandline(Input.argparse_settings())["config_file"]

# Set up and run the coupled simulation
cs = CoupledSimulation(config_file)
run!(cs)

# Postprocessing
if !skip_postprocess
    # Load the Makie plotting stack now (unless it was already loaded above) so the
    # simulation reaches its first step without compiling the plotting packages.
    plots_during_run || include(joinpath(@__DIR__, "..", "load_plotting.jl"))
    conservation_softfail =
        Input.get_coupler_config_dict(config_file)["conservation_softfail"]
    rmse_check = Input.get_coupler_config_dict(config_file)["rmse_check"]
    postprocess(cs; conservation_softfail, rmse_check)
end
