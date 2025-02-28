# ClimaCoupler.jl

ClimaCoupler.jl provides coupled system time stepping control and support for mapping import and export
boundary information between components. It can handle atmosphere, ocean, land, and sea ice component models,
and the source code is agnostic to the internals of these component models.

ClimaCoupler.jl contains a directory `experiments/ClimaEarth/` which contains
additional infrastructure to run coupled AMIP experiments using component models
from ClimaAtmos.jl and ClimaLand.jl. This includes machinery for visualizing output
in the `user_io/` folder, and component model-specific initialization, access,
and other helper functions in the `components/` folder, which will soon be moved to
the respective component model packages.

Additional smaller coupling examples can be found in the `experiments/ClimaCore/` directory.
These are meant to serve as an introduction to coupling and the types of functionality
required for it.

<!-- Links and shortcuts -->
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/ClimaCoupler.jl/dev/

[docs-bld-img]: https://github.com/CliMA/ClimaCoupler.jl/workflows/Documentation/badge.svg
[docs-bld-url]: https://github.com/CliMA/ClimaCoupler.jl/actions?query=workflow%3ADocumentation

[unit-tests-img]: https://github.com/CliMA/ClimaCoupler.jl/actions/workflows/ci.yml/badge.svg
[unit-tests-url]: https://github.com/CliMA/ClimaCoupler.jl/actions?query=workflow%3Aci

[codecov-img]: https://codecov.io/gh/CliMA/ClimaCoupler.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/ClimaCoupler.jl

[downloads-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FClimaCoupler&query=total_requests&suffix=%2Ftotal&label=Downloads
[downloads-url]: http://juliapkgstats.com/pkg/ClimaCoupler

|||
|---------------------:|:-----------------------------------------------|
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url] [![docs build][docs-bld-img]][docs-bld-url]|
| **Unit Tests**       | [![unit tests][unit-tests-img]][unit-tests-url] [![codecov][codecov-img]][codecov-url]|
| **Downloads**        | [![downloads][downloads-img]][downloads-url]|
|||

Recommended Julia Version: Stable release v1.11.1. CI tests Julia v1.10 and 1.11.

# Running AMIP
Here we will focus on the AMIP experiment, which uses the environment in the `experiments/ClimaEarth/` subdirectory of ClimaCoupler.jl
The first step to do this is to install all required packages for the environment using the following Julia command:
```julia
julia --project=experiments/ClimaEarth -E "using Pkg; Pkg.instantiate(); Pkg.build()"
```

Now you're ready to run the experiment, which uses the `run_amip.jl` driver. To run interactively:
```
julia --project=experiments/ClimaEarth
julia> include("experiments/ClimaEarth/run_amip.jl")
```

Or to run from the terminal:
```julia
julia --project=experiments/ClimaEarth experiments/ClimaEarth/run_amip.jl
```

When running from the terminal, you can also specify a configuration file to use for the simulation setup, and a job ID to keep track of this run's output.
Existing configuration files are specified in the `config/` directory within ClimaCoupler.jl.
For example, to run the default AMIP configuration, you could use the following command:
```julia
julia --project=experiments/ClimaEarth experiments/ClimaEarth/run_amip.jl --config_file config/ci_configs/amip_default.yml --job_id amip_default
```

Output from your run will be saved in the folder `experiments/ClimaEarth/output/<job_id>/`.
If no configuration file is specified, the default
`amip_default.yml` will be used, and output will be saved in `experiments/ClimaEarth/output/amip_default/`.

The output will take up approximately 1GB of space, and the simulation will take around 10 minutes to run on a single CPU, or less time on multiple CPUs or GPU.

Note: If you want to set the configuration file to something other than the default
while running the driver interactively, you'll need to
manually set the value for `config_file`.

For example, to use the configuration file found at `config/ci_configs/amip_default.yml`, you would set `config_file` as follows in the `run_amip` driver:
```
config_file = "config/ci_configs/amip_default.yml"
```

### A Note about ClimaComms and MPI
If you don't intend to run your simulation using MPI, but you see an error about MPI and your simulation crashes,
ClimaComms may be incorrectly selecting the configuration for your run.
In this case, you can force ClimaComms to ignore MPI with
```
export CLIMACOMMS_CONTEXT="SINGLETON"
```
from the terminal, or
```
ENV["CLIMACOMMS_CONTEXT"]="SINGLETON"
```
from within the Julia environment before running the experiment.

Sometimes this happens when you are running in an interactive SLURM session.

## Running on GPU or with MPI

### CUDA.jl and MPI.jl packages
CUDA.jl and MPI.jl are required to run ClimaCoupler's AMIP experiment on GPU and in parallel (with MPI), respectively.
Not every machine is capable of running on GPU or with MPI, and the AMIP experiment can be run on CPU
without these packages, so they aren't included in the ClimaCoupler AMIP experiment environment.
This means that if you want to run our driver using these capabilities (locally or remotely), you should install
CUDA and MPI in your machine's base Julia environment. This will make the packages available to the
AMIP experiment environment, as it is a sub-environment of the base environment.
If you wish to run the AMIP experiment with GPU or MPI locally as well as remotely,
this process must be done on each machine you want to run on.

This can be done using the following command in the terminal from any directory:
```julia
julia -E "using Pkg; Pkg.add(\"CUDA\"); Pkg.add(\"MPI\")"
```

Now, if you enter your base environment by running `julia` and then check the packages with `] st`, you should see something like:
```
(@v1.10) pkg> st
Status `~/.julia/environments/v1.10/Project.toml`
  [<hash>] CUDA vX.Y.Z
  [<hash>] MPI vX.Y.Z
```

### Environment variables
Additionally, there are some environment variables we must set in these cases.

To run on GPU, we need to run `export CLIMACOMMS_DEVICE="CUDA"` in the terminal, or
`ENV["CLIMACOMMS_DEVICE]="CUDA"` within the Julia environment _before_ running the experiment.

To run with MPI, we need to run `export CLIMACOMMS_CONTEXT="MPI"` in the terminal, or
`ENV["CLIMACOMMS_CONTEXT]="MPI"` within the Julia environment _before_ running the experiment.

## Caltech users: Running AMIP remotely
The main difference between running code locally vs running remotely is
the module loading step. CliMA uses [ClimaModules](https://github.com/CliMA/ClimaModules?tab=readme-ov-file#clima-modules-for-new-central) to coordinate the modules
needed to run CliMA code on Caltech's clusters.

On Central, you can load the appropriate module package by running the following in the terminal:
```
export MODULEPATH="/groups/esm/modules:$MODULEPATH"
module load climacommon
```

> Remember: This should be done _after_ requesting a compute node, using the command `srun --pty -t 01:00:00 -p expansion bash` or similar

On `clima`, you can load the appropriate module package by running the following in the terminal:
```
module load common
```

For additional information about these clusters, including how to gain access for the first time,
see our slurm-buildkite wiki pages for [Central](https://github.com/CliMA/slurm-buildkite/wiki/Central) and [clima](https://github.com/CliMA/slurm-buildkite/wiki/clima).

# Running Slabplanet
The `run_amip.jl` driver contains two modes: the full AMIP mode and a Slabplanet mode, where all surfaces are thermal slabs. Since AMIP is not a closed system, the Slabplanet mode is useful for checking conservation properties of the coupling.

Running a Slabplanet simulation is the same as running an AMIP simulation, except for the specifics of the configuration file provided, so all information from the `Running AMIP` section will apply here too. Note that the default configuration used by `run_amip.jl` specifies an AMIP simulation, so a configuration file must be specified to run a Slabplanet simulation. This can be done as follows:
```julia
julia --project=experiments/ClimaEarth experiments/ClimaEarth/run_amip.jl --config_file config/ci_configs/slabplanet_default.yml --job_id slabplanet_default
```

To ensure that conservation is tracked throughout the experiment, the `energy_check` field of the configuration file must be set to true.
