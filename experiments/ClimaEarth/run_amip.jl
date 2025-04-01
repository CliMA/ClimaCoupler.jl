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
julia --project=experiments/ClimaEarth experiments/ClimaEarth/run_amip.jl
```

## Configuration
You can also specify a custom configuration file to run the coupled simulation
in a different setup. The configuration file should be a TOML file that overwrites
the input fields specified in `experiments/ClimaEarth/cli_options.jl`.
A set of example configuration files can be found in the `config/ci_configs/` directory.

For example, to run the coupled simulation with a different configuration file:
```bash
julia --project=experiments/ClimaEarth experiments/ClimaEarth/run_amip.jl --config_file="path/to/config.toml"
```

To run the coupled simulation interactively with a different configuration file,
set the `config_file` variable in this script to be the path to that file.

For more details about running a coupled simulation, including how to run in a
Slabplanet configuration, please see our [README.md](https://github.com/CliMA/ClimaCoupler.jl/blob/main/README.md).
=#

# Load the necessary modules and code to run the coupled simulation
include("setup_run.jl")

# Get the configuration file from the command line (or manually set it here)
config_file = parse_commandline(argparse_settings())["config_file"]

# Set up and run the coupled simulation
cs = CoupledSimulation(config_file)
run!(cs)

# Postprocessing
# TODO: Remove this option?
conservation_softfail = get_coupler_config_dict(config_file)["conservation_softfail"]
postprocess(cs, conservation_softfail)
