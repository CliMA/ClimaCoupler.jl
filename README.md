<div align="center">
  <img src="docs/src/assets/logo.svg" alt="ClimaCoupler.jl Logo" width="128" height="128">
</div>

# ClimaCoupler.jl

ClimaCoupler.jl provides coupled system time stepping control and support for mapping import and export
boundary information between components. It can handle atmosphere, ocean, land, and sea ice component models,
and the source code is agnostic to the internals of these component models.

ClimaCoupler.jl targets global, long-term climate simulations, specifically for the
[Coupled Model Intercomparison Project (CMIP)](https://wcrp-cmip.org) and
[Atmospheric Model Intercomparison Project (AMIP)](https://pcmdi.llnl.gov/mips/amip/home/overview.html).
With these objectives in mind, this coupler is also designed to be flexible and modular,
so it supports simpler slabplanet-type setups as well.

|                    |                                                                                    |                    |                                         |
|-------------------:|:-----------------------------------------------------------------------------------|-------------------:|:----------------------------------------|
| **Documentation** | [![stable][docs-stable-img]][docs-stable-url] [![dev][docs-dev-img]][docs-dev-url] | **License**        | [![license][license-img]][license-url]  |
| **Tests**         | [![gha ci][gha-ci-img]][gha-ci-url] [![buildkite][bk-ci-img]][bk-ci-url]           | **Code Coverage**  | [![codecov][codecov-img]][codecov-url]  |
| **Downloads**     | [![Downloads][dlt-img]][dlt-url]                                                   | **Nightly AMIP**   | [![buildkite][bk-amip-img]][bk-amip-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://CliMA.github.io/ClimaCoupler.jl/stable/

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/ClimaCoupler.jl/dev/

[license-img]: https://img.shields.io/badge/license-Apache%202.0-blue.svg
[license-url]: https://github.com/CliMA/ClimaCoupler.jl/blob/main/LICENSE

[gha-ci-img]: https://github.com/CliMA/ClimaCoupler.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/ClimaCoupler.jl/actions/workflows/ci.yml

[bk-ci-img]: https://badge.buildkite.com/1d6cfcc4219656239b9917210fa4c6e6f82b0516c2b0074fe1.svg?branch=main
[bk-ci-url]: https://buildkite.com/clima/climacoupler-ci

[bk-amip-img]: https://badge.buildkite.com/67ed8ae5c42ef2ea566fcf3a981e0165b39600dd24c5f9a37d.svg?branch=main
[bk-amip-url]: https://buildkite.com/clima/climacoupler-coarse-nightly-amip

[codecov-img]: https://codecov.io/gh/CliMA/ClimaCoupler.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/ClimaCoupler.jl

[dlt-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FClimaCoupler&query=total_requests&suffix=%2Ftotal&label=Downloads
[dlt-url]: https://juliapkgstats.com/pkg/ClimaCoupler

## Running AMIP

### Setting up

Clone the repository and run every command below from its root directory. Config file paths
and the `output/` directory are resolved relative to the directory Julia is started from.
```bash
git clone https://github.com/CliMA/ClimaCoupler.jl.git
cd ClimaCoupler.jl
```

Use Julia 1.10 or 1.11. The experiment environments ship manifests pinned to those versions
(`Manifest.toml` for 1.10, `Manifest-v1.11.toml` for 1.11); other Julia versions fall back to
the 1.10 manifest and may need to re-resolve it. With [juliaup](https://github.com/JuliaLang/juliaup),
`juliaup add 1.11` installs a matching version and `julia +1.11` selects it.

AMIP, slabplanet, and subseasonal runs use the environment in `experiments/AMIP/`. Install its
packages once, and again after pulling changes to its manifest:
```julia
julia --project=experiments/AMIP -E "using Pkg; Pkg.instantiate(); Pkg.build()"
```

The first simulation also downloads its input data (land-sea masks, SST and sea-ice data, etc.)
as Julia artifacts.

### Running from the terminal

```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl
```

This runs `config/ci_configs/amip_default.yml`, a short test configuration, and saves output to
`output/amip_default/`. To choose a configuration and a job ID for the output folder:
```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl --config_file config/ci_configs/amip_default.yml --job_id my_amip_run
```

Output is saved in `output/<job_id>/`; the job ID defaults to the config file name.
Ready-to-use configuration files are in `config/ci_configs/`, and all configuration options
are listed on the [Input](https://CliMA.github.io/ClimaCoupler.jl/dev/input/) docs page.
Most of the run time of a first launch is compilation.

### Running from the REPL

Start Julia from the repository root with `julia --project=experiments/AMIP`, then load the
packages and pass a configuration file to `CoupledSimulation`:
```julia
include("experiments/AMIP/code_loading.jl")

cs = CoupledSimulation("config/ci_configs/amip_default.yml")
run!(cs)
postprocess(cs)
```

To change the setup, write your own YAML file in `config/`, or load an existing one into a
dictionary, edit it, and pass the dictionary instead:
```julia
config_dict = Input.get_coupler_config_dict("config/ci_configs/amip_default.yml")
config_dict["t_end"] = "1days"

cs = CoupledSimulation(config_dict)
run!(cs)
postprocess(cs)
```

The [Running a Simulation](https://CliMA.github.io/ClimaCoupler.jl/dev/running/) docs page also
covers stepping a simulation manually and using ClimaCoupler as an installed package.

### Running CMIP

CMIP runs use their own environment, `experiments/CMIP/`, which adds Oceananigans and
ClimaSeaIce, and must be given a CMIP configuration file:
```bash
julia --project=experiments/CMIP -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/CMIP experiments/CMIP/run_simulation.jl --config_file config/ci_configs/cmip_oceananigans_climaseaice.yml
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

### Running on GPU or with MPI

Both experiment environments already include CUDA.jl and MPI.jl; `code_loading.jl` loads them
through `ClimaComms.@import_required_backends` when the environment variables below request them.
Set the variables before starting Julia, or with `ENV` _before_ including `code_loading.jl`:
```bash
export CLIMACOMMS_DEVICE="CUDA"   # run on an NVIDIA GPU
export CLIMACOMMS_CONTEXT="MPI"   # run distributed, with one Julia process per MPI rank
```
or, in the REPL:
```julia
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "MPI"
include("experiments/AMIP/code_loading.jl")
```

Keep the `device` configuration option at its default, `"auto"`, so that `CLIMACOMMS_DEVICE`
selects the device. For CUDA driver compatibility and checking that a run is really on the GPU,
see [running_on_gpu.md](docs/dev-guides/workflow/running_on_gpu.md).

### Caltech users: Running AMIP remotely
The main difference between running code locally vs running remotely is
the module loading step. CliMA uses [ClimaModules](https://github.com/CliMA/ClimaModules?tab=readme-ov-file#clima-modules-for-new-central) to coordinate the modules
needed to run CliMA code on Caltech's clusters.

On Central, you can load the appropriate module package by running the following in the terminal:
```
export MODULEPATH="/resnick/groups/esm/modules:$MODULEPATH"
module load climacommon
```

> Remember: This should be done _after_ requesting a compute node, using the command `srun --pty -t 01:00:00 -p expansion bash` or similar

On `clima`, you can load the appropriate module package by running the following in the terminal:
```
module load common
```

For additional information about these clusters, including how to gain access for the first time,
see our slurm-buildkite wiki pages for [Central](https://github.com/CliMA/slurm-buildkite/wiki/Central) and [clima](https://github.com/CliMA/slurm-buildkite/wiki/clima).

## Running Slabplanet
The `run_simulation.jl` driver contains two modes: the full AMIP mode and a Slabplanet mode, where all surfaces are thermal slabs. Since AMIP is not a closed system, the Slabplanet mode is useful for checking conservation properties of the coupling.

Running a Slabplanet simulation is the same as running an AMIP simulation, except for the specifics of the configuration file provided, so all information from the [Running AMIP](#running-amip) section will apply here too. Note that the default configuration used by `run_simulation.jl` specifies an AMIP simulation, so a configuration file must be specified to run a Slabplanet simulation. This can be done as follows:
```julia
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl --config_file config/ci_configs/slabplanet_default.yml --job_id slabplanet_default
```

To ensure that conservation is tracked throughout the experiment, the `energy_check` field of the configuration file must be set to true.

## Contributing

Contributors should follow the shared CliMA engineering standards in [`docs/dev-guides/`](docs/dev-guides/), which cover architecture, performance, code quality, documentation, and workflows. These are vendored from [CliMA/DeveloperGuides](https://github.com/CliMA/DeveloperGuides). The repo's [`AGENTS.md`](AGENTS.md) is a starting point for AI agents with repo-specific guidance.
