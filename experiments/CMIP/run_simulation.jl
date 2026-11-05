# # CMIP Driver

#=
## Overview

CMIP is a standard experimental protocol of the World Climate Research Programme (WCRP).
It is used as a model benchmark for coupled, prognostic atmospheric, land, ocean, and sea ice model components.

For more information, see the WCRP's specifications for [CMIP](https://wcrp-cmip.org).

## Running the CMIP configuration
To run a coupled simulation in the default CMIP configuration, run the
following command from the root directory of the repository:
```bash
julia --project=experiments/CMIP experiments/CMIP/run_simulation.jl --config_file config/ci_configs/cmip_oceananigans_climaseaice.yml
```

## Configuration
You can also specify a custom configuration file to run the coupled simulation
in a different setup. The configuration file is a YAML file that overwrites
the input fields specified in the ClimaCoupler Input module.
A set of example configuration files can be found in the `config/ci_configs/` directory.

To run the coupled simulation interactively, start Julia from the root directory with
`julia --project=experiments/CMIP`, run `include("experiments/CMIP/code_loading.jl")`,
and pass a configuration file path to `CoupledSimulation`. To change settings, edit a YAML
file, or load one with `Input.get_coupler_config_dict`, modify the returned dictionary,
and pass the dictionary to `CoupledSimulation`.

For more details about running a coupled simulation, including how to run in a
Slabplanet configuration, please see our documentation.
=#

# Load the necessary modules to run the coupled simulation
include("code_loading.jl")

# Get the configuration file from the command line
config_file = Input.parse_commandline(Input.argparse_settings())["config_file"]

# Set up and run the coupled simulation
cs = CoupledSimulation(config_file)
run!(cs)

# Postprocessing
conservation_softfail = Input.get_coupler_config_dict(config_file)["conservation_softfail"]
rmse_check = Input.get_coupler_config_dict(config_file)["rmse_check"]
postprocess(cs; conservation_softfail, rmse_check)
