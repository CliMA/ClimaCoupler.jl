# Simulations & Timestepping

## Simulations Interface
`ClimaCoupler.jl` organizes coupled models and their execution
through the use of `Simulation`s. A `Simulation` wraps a physical model and
its time integrator; it contains all necessary information to execute a 
model run. Component models subtype [`ClimaCoupler.AbstractSimulation`](@ref)
so that coupling methods can be implmented and dispatch by model type.
Component `Simulation`s are collected in a [`CoupledSimulation`](@ref) along with
the coupler itself. This is the full coupled modeling system. 

When coupling, component models must implement three methods:
- [`step!`](@ref): advances the component model by a specified step size
- [`coupler_push!`](@ref): prepares and puts coupled fields from the model
    into the coupler via `coupler_put!` calls
- [`coupler_pull!`](@ref): gets coupled fields from the coupler
    via `coupler_get!` calls and prepares them to be ingested by the model.
These methods hide each components backend implementations (for timestepping or 
field access) from the coupler, allowing heterogeneity in model backends and
the standardization of some coupled timestepping schemes. However, like any
other simulation, a `CoupledSimulation` may also implement a
[`ClimaCoupler.step!`](@ref) method, which in this context specifies
the coupling scheme details (e.g. explicit, leap-frog, concurrent, etc...).

```@docs
ClimaCoupler.AbstractSimulation
ClimaCoupler.CoupledSimulation
ClimaCoupler.run!
ClimaCoupler.step!
```