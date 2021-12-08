# Coupled Simulations & Timestepping

`ClimaCoupler` organizes coupled models and their execution
via the [`CoupledSimulation`](@ref) interface. An implementation
of a [`CoupledSimulation`](@ref), like any CliMA simulation, must 
implement a [`ClimaCoupler.step!`](@ref) method, which in this context specifies
the coupling scheme details (e.g. explicit, leap-frog, concurrent, etc...).

Component models being coupled must provide three methods for use in a 
[`CoupledSimulation`](@ref)'s [`ClimaCoupler.step!`](@ref) method:
- `step!`: advances the component model in time
- [`coupler_push!`](@ref): prepares and puts coupled fields from the model
    into the coupler via `coupler_put!` calls
- [`coupler_pull!`](@ref): gets coupled fields from the coupler
    via `coupler_get!` calls and prepares them to be ingested by the model.

```@docs
ClimaCoupler.CoupledSimulation
ClimaCoupler.run!
ClimaCoupler.step!
ClimaCoupler.coupler_push!
ClimaCoupler.coupler_pull!
```