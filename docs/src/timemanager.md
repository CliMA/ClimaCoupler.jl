# TimeManager

This module contains functions that handle dates and times
in simulations. The functions in this module often call
functions from Julia's [Dates](https://docs.julialang.org/en/v1/stdlib/Dates/) module.

## TimeManager API

```@docs
TimeManager.maybe_trigger_callback
```

## ITime

`ITime`, or _integer time_, is a time type used by CliMA simulations to keep
track of simulation time. For more information, refer to the
[TimeManager section](https://clima.github.io/ClimaUtilities.jl/dev/timemanager/)
in ClimaUtilities and the [ITime section](https://clima.github.io/ClimaAtmos.jl/dev/itime/)
in ClimaAtmos.

### How do I use ITime?

If you are running a simulation from a YAML file, you can simply set `use_itime`
to true to enable `ITime`. If you do not want to use `ITime` and want to use
floating point numbers, then set `use_itime` to false to not use `ITime`.
