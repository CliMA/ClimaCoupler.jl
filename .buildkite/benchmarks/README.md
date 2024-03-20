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

#### CPU ClimaAtmos with diagnostic EDMF
- Atmosphere-only simulation
- Run on 64 CPU threads

#### CPU AMIP with diagnostic EDMF
- ClimaAtmos coupled to ClimaLand bucket model, with prescribed sea surface
temperature and sea ice
- Run on 64 CPU threads

#### GPU ClimaAtmos with diagnostic EDMF
- Atmosphere-only simulation
- Run on 4 A100 GPUs sharing 1 node

#### GPU AMIP with diagnostic EDMF
- ClimaAtmos coupled to ClimaLand bucket model, with prescribed sea surface
temperature and sea ice
- Run on 4 A100 GPUs sharing 1 node

### Comparison Metrics
- Simulated years per day (SYPD): The number of years of simulation time we
can run in 1 day of walltime
- CPU simulation object allocations: The allocations in GB of the simulation
object, which contains everything needed to run the simulation.
In the atmosphere-only case, this is the `AtmosSimulation` object.
In the coupled case, this is the `CoupledSimulation` object, which includes
all of the component models, coupler fields, and auxiliary objects. More
information on this object can be found in the `Interfacer` docs.
