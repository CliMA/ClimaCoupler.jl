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
in a different setup. The configuration file should be a TOML file that overwrites
the input fields specified in the ClimaCoupler Input module.
A set of example configuration files can be found in the `config/ci_configs/` directory.

To run the coupled simulation interactively with a different configuration file,
set the `config_file` variable in this script to be the path to that file.

For more details about running a coupled simulation, including how to run in a
Slabplanet configuration, please see our documentation.
=#

# Load the necessary modules to run the coupled simulation
include("code_loading.jl")

# Get the configuration file from the command line (or manually set it here)
config_file = Input.parse_commandline(Input.argparse_settings())["config_file"]

# Set up and run the coupled simulation
cs = CoupledSimulation(config_file)


# TODO: Remove this before merging

import Oceananigans as OC

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model

    (Tmax, Tmin) = extrema(ocean.tracers.T)
    (Smax, Smin) = extrema(ocean.tracers.S)
    (ηmax, ηmin) = extrema(ocean.free_surface.displacement)
    umax = maximum(ocean.velocities.u)
    vmax = maximum(ocean.velocities.v)
    wmax = maximum(ocean.velocities.w)
    step_time = 1e-9 * (time_ns() - wall_time[])
    @info "time: $(OC.prettytime(sim)), iteration: $(OC.iteration(sim)), Δt: $(OC.prettytime(sim.Δt)), " *
          "extrema(η): ($(round(ηmin, sigdigits=2)), $(round(ηmax, sigdigits=2))) " *
          "extrema(T, S): ($(round(Tmin, digits=2)), $(round(Tmax, digits=2))) ᵒC, " *
          "($(round(Smin, digits=2)), $(round(Smax, digits=2))) psu " *
          "maximum(u): ($(round(umax, sigdigits=2)), $(round(vmax, sigdigits=2)), $(round(wmax, sigdigits=2))) m/s, " *
          "wall time: $(OC.prettytime(step_time))"

    wall_time[] = time_ns()

    return nothing
end

# Attaching a progress function to the ocean
OC.add_callback!(cs.model_sims.ocean_sim.ocean, progress, OC.IterationInterval(1))

run!(cs)

# Postprocessing
conservation_softfail = Input.get_coupler_config_dict(config_file)["conservation_softfail"]
rmse_check = Input.get_coupler_config_dict(config_file)["rmse_check"]
postprocess(cs; conservation_softfail, rmse_check)
