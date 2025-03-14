## ClimaCoupler Benchmarks Pipeline

### Purpose
The goal of the benchmarks pipeline is to have concrete comparisons between
analogous simulations of different setups and on different architectures.
This allows us to compare things like performance and allocations across
atmosphere-only vs coupled runs, and on CPU vs GPU.

This pipeline is triggered manually rather than on a schedule, so that we
can monitor the various metrics after specific changes made to the code.

### Simulation Setups
#### All simulations
- Timestep: 120 seconds
- Horizontal resolution: 30 spectral elements (~110km)
- Vertical resolution: 63 levels
- Config setup duplicated from ClimaAtmos.jl v0.23.0
[gpu_aquaplanet_diagedmf.yml](https://github.com/CliMA/ClimaAtmos.jl/blob/v0.23.0/config/gpu_configs/gpu_aquaplanet_diagedmf.yml),
with minor tweaks

#### CPU ClimaAtmos without diagnostic EDMF
- Atmosphere-only simulation
- Run on 64 CPU threads
- Frierson diffusion model

#### CPU ClimaAtmos with diagnostic EDMF
- Atmosphere-only simulation
- Run on 64 CPU threads
- Diagnostic EDMF parameterization

#### CPU AMIP with diagnostic EDMF
- ClimaAtmos coupled to ClimaLand bucket model, with prescribed sea surface
temperature and sea ice
- Run on 64 CPU threads

#### CPU AMIP with diagnostic EDMF and IO
- ClimaAtmos coupled to ClimaLand bucket model, with prescribed sea surface
temperature and sea ice
- Run on 64 CPU threads
- Includes diagnostics output every 10 hours

#### GPU ClimaAtmos without diagnostic EDMF
- Atmosphere-only simulation
- Run on 4 A100 GPUs sharing 1 node
- Frierson diffusion model

#### GPU ClimaAtmos with diagnostic EDMF
- Atmosphere-only simulation
- Run on 4 A100 GPUs sharing 1 node
- Diagnostic EDMF parameterization

#### GPU AMIP with diagnostic EDMF
- ClimaAtmos coupled to ClimaLand bucket model, with prescribed sea surface
temperature and sea ice
- Run on 4 A100 GPUs sharing 1 node

#### GPU AMIP with diagnostic EDMF and IO
- ClimaAtmos coupled to ClimaLand bucket model, with prescribed sea surface
temperature and sea ice
- Run on 4 A100 GPUs sharing 1 node
- Includes diagnostics output every 10 hours

### Comparison Metrics
- Simulated years per day (SYPD): The number of years of simulation time we
can run in 1 day of walltime
- CPU maximum resident set size (max RSS): The max RSS memory footprint on the
CPU of this process since it began. This is measured for both CPU and GPU runs.
