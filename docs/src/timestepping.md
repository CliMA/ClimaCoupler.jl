# Coupled Timestepping

`CouplerMachine` currently uses a sequential timestepping framework in which one
component steps forward before passing its updated state to another. An outer
"coupled" timestep determines when component models synchronize and coordinate
with the coupler. Within this coupled timestep, components take an integer number
of substeps, and evolve independently from each other.

`CplSolver` extends the ODE solver API of
[ClimateMachine.ODESolvers](https://clima.github.io/ClimateMachine.jl/latest/APIs/Numerics/ODESolvers/ODESolvers/).

```@docs
CouplerMachine.CplSolver
```