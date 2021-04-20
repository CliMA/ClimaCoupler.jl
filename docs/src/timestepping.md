# Coupled Timestepping

`CouplerMachine` uses a concurrent timestepping framework. An outer "coupled" timestep
determines when component models synchronize and coordinate with the coupler.
Within this coupled timestep, components take an integer number of substeps, and evolve
independently from each other.

| ![Coupled Timestepping](images/cpltimestep.png) |
|:--:|
| *Coupled timestepping with two component models.* |

`CplSolver` extends the ODE solver API of
[ClimateMachine.ODESolvers](https://clima.github.io/ClimateMachine.jl/latest/APIs/Numerics/ODESolvers/ODESolvers/).

```@docs
CouplerMachine.CplSolver
```