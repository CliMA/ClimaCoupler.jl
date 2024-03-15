# ClimaCoupler.jl

Coupler Specific Shared Development

<!-- Links and shortcuts -->
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/ClimaCoupler.jl/dev/

[docs-bld-img]: https://github.com/CliMA/ClimaCoupler.jl/workflows/Documentation/badge.svg
[docs-bld-url]: https://github.com/CliMA/ClimaCoupler.jl/actions?query=workflow%3ADocumentation

[unit-tests-img]: https://github.com/CliMA/ClimaCoupler.jl/workflows/Unit%20Tests/badge.svg
[unit-tests-url]: https://github.com/CliMA/ClimaCoupler.jl/actions?query=workflow%3A%22Unit+Tests%22

[codecov-img]: https://codecov.io/gh/CliMA/ClimaCoupler.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/ClimaCoupler.jl

|||
|---------------------:|:-----------------------------------------------|
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url] [![docs build][docs-bld-img]][docs-bld-url]|
| **Unit Tests**       | [![unit tests][unit-tests-img]][unit-tests-url] [![codecov][codecov-img]][codecov-url]|

Provides coupled system time stepping control and support for mapping import and export
boundary information between components.

Recommended Julia Version: Stable release v1.10.1. CI no longer tests earlier versions of Julia.

## Start Up
Before starting Julia, ensure your environment is properly set up. This includes loading the correct Julia version (specified in the `Project.toml` file), modules and setting the correct environment variables.

Typically (e.g., when running on your local machine) the installation of the required module libraries and specification of environment variables should be automatically handled by the Julia package manager upon running `]instantiate` and `]build` in the Julia REPL,
which reads the `Project.toml` file. Each experiment has its own `Project.toml` file and `Manifest.toml` file, which are used to specify the exact environment for the experiment. Once instantiated, the experiment run scripts should be able to run without any further environment setup (e.g., with `julia --project --threads 8 experiments/<your_experiment>.jl`).

For computationally intensive jobs, it is recommended to run the jobs on a cluster using MPI and/or GPU triggered by a SLURM job script.
In this case, we need to load specific MPI and GPU compatible modules and set the correct environment variables. For an example of a SLURM job script, see `test/mpi_tests/local_checks.sh`.

