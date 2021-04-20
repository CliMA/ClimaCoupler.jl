# Coupler Object

The CouplerMachine defines a type [`CplState`](@ref) for a _container_ variable
that holds information about the field boundary values that are being used to
couple components. Components can use a [`CouplerMachine.put!`](@ref) operation to 
export a set of field values to a `CplState` variable. A [`CouplerMachine.get`](@ref)
operation is used to retrieve a set field values from a `CplState` variable.
During this exchange, the coupler manages ancillary operations such as 
regridding, unit conversions, or filtering.

## Coupler Object API

```@docs
    CouplerMachine.CplState
    CouplerMachine.register_cpl_field!
    CouplerMachine.put!
    CouplerMachine.get
```
