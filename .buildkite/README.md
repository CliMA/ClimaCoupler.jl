# ClimaCoupler.jl Buildkite Pipelines

ClimaCoupler.jl has a number of Buildkite pipelines that are run regularly.
More information about each of them is provided here.


# [ClimaCoupler-CI](https://buildkite.com/clima/climacoupler-ci)
- Triggered for each PR; must pass for PR to be merged
- Walltime (approx.): 45 minutes
- Runs on Caltech Central cluster
- Configuration file location: `config/ci_configs/`

This pipeline focuses on running unit tests in different setups, as well as small
integration tests with different configurations. This includes running the
ClimaCoupler.jl test suite on GPU,
some unit tests with MPI, and short coupled AMIP and slabplanet simulations.

Most of the coupled simulations run for a few timesteps (300-800 seconds), but
there are some longer ones that run for 10 days, and a coarse run that runs for
30 days of simulation time. The resolution of these simulations ranges from 4 to
16 horizontal elements. The simulations test all options for turbulent flux partitioning,
float types, albedo specification, and land surface types. Simulations are run on both
CPU and GPU; those run on GPU are labeled as such.

We try to limit runtime of simulations in this pipeline to minimize its overall runtime
and increase our development speed. As a general rule, it's best not to add any run
that increases the overall runtime of this pipeline. This can most easily be achieved by
reducing simulation length or resolution.

## [ClimaCoupler-ClearDepot](https://buildkite.com/clima/climacoupler-cleardepot)
The `ClimaCoupler-CI` pipeline uses a shared depot to cache some packages and reduce precompilation time.
The depot is a folder that contains downloaded and compiled julia packages and artifacts and it is shared across
different runs of the same pipeline. While this is helpful to decrease walltime of that pipeline, the depot sometimes
gets corrupted and needs to be reset. The `ClimaCoupler-ClearDepot` pipeline can be
manually triggered to reset the depot.

# [ClimaCoupler - Coarse Nightly AMIP](https://buildkite.com/clima/climacoupler-coarse-nightly-amip)
- Scheduled daily Sunday-Thursday at 6pm PST (2am UTC +1 day)
- Walltime (approx.): 12 hours
- Runs on Caltech Clima node (4 GPUs)
- Configuration file location: `config/nightly_configs/`

This pipeline runs nightly on weekdays and uses one configuration which specifies
a coarse AMIP running for 1+ years. The pipeline triggers 4 runs using this configuration,
each on 1 GPU.
Each run sets its own random seed so the simulations vary slightly, which allows us to more
robustly test the stability of the simulations. The random seed of each simulation
is printed so they can be reproduced, and 1 of the 4 runs always uses the same seed.

This pipeline runs with the main branches of all major CliMA packages and its main goal is
the goal to monitor stability of a coarse AMIP simulation and to catch potential problems
early. This pipeline is scheduled, but it can also be used to test out configuration changes
in a smaller test case before running higher resolution global AMIP runs.

> _Note:_  This configuration is under development and frequently changing.

# [ClimaCoupler - CPU/GPU Benchmarks](https://buildkite.com/clima/climacoupler-cpu-gpu-benchmarks)
- Scheduled weekly Sunday at 12am PST (8am UTC)
- Walltime (approx.): 16 hours
- Runs on Caltech Clima node (all GPUs)
- Configuration file location: `config/benchmarks_configs/`

This pipeline includes the following runs, each on both CPU and GPU:
- ClimaAtmos without diagnostic EDMF
- ClimaAtmos with diagnostic EDMF
- Coupled AMIP with diagnostic EDMF
- Coupled AMIP with diagnostic EDMF and IO

Each of the runs has a simulation length of 12 hours, and a resolution of 30
elements in the horizontal and 63 in the vertical with a 55km atmosphere top.

The goal of this pipeline is to be able to compare performance of uncoupled atmosphere
runs vs coupled AMIP runs, and to see the performance impact of slow things like
computationally-intensive parameterizations or diagnostics.

This pipeline also produces a table showing useful comparison metrics for each of
the above runs, such as SYPD or maximum memory usage. This table is sent to the
CliMA Slack `#coupler-report` channel upon successful completion of all runs.

For more information, see [.buildkite/benchmarks/README.md](benchmarks/README.md).

# [ClimaCoupler-AMIP](https://buildkite.com/clima/climacoupler-amip)
- Scheduled weekly Monday at 8pm PST (4am UTC +1 day)
- Walltime (approx.): 80 hours
- Runs on Caltech Clima node (1 GPU)
- Configuration file location: `config/amip_configs/`

This pipeline runs the target AMIP configuration that we're currently working to stabilize.
The simulation runs for 3 years and has a resolution of 16
elements in the horizontal and 63 in the vertical
and includes diagnostic EDMF in the atmosphere. It is run on 1 GPU.

# [ClimaCoupler-LongRuns](https://buildkite.com/clima/climacoupler-longruns)
- Scheduled weekly Sunday at 12am PST (8am UTC)
- Walltime (approx.): 24 hours
- Runs on Caltech Central cluster
- Configuration file location: `config/longruns/`

This pipeline provides a set of simulations of increasing physical and computational complexity.

The simplest setup is a `slabplanet_aqua` simulation, running an atmosphere model with a constant-temperature
thermal slab ocean model as the only surface model, and
with each model running independently
(i.e. no flux exchange by the coupler). Complexity increases include
evolving the surface model, adding in flux exchange, using multiple surface models,
using more complex parameterizations (e.g. diagnostic EDMF), and increasing spatial resolution.
Including these degrees of complexity one at a time allows us to track the stability of our
coupled model across different setups.

The resolutions of runs in this pipeline vary from 4 to 16 horizontal elements.
The simplest simulations are only run on CPU, but some of the more complex
simulations are run on GPU in addition to CPU.

