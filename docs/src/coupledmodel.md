# Coupled Model Components

CouplerMachine provides a wrapper for model components so that they may
connect to the coupler. The [`CplModel`](@ref) struct packages a component
model with the information needed for the [`CplSolver`](@ref) to evolve it.

```@docs
    CouplerMachine.CplModel
```