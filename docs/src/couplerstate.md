# Coupler State

The coupler provides a space to store information being passed between coupled model components at their boundaries. During this exchange, the coupler manages ancillary operations such as regridding, unit conversions, filtering, etc.

The ClimaCoupler defines a type [`CouplerState`](@ref) for a _container_ variable
that holds information about the field boundary values that are being used to
couple components. Components can use a [`coupler_put!`](@ref) operation to 
export a set of field values to a `CouplerState` variable. A [`coupler_get`](@ref)
operation is used to retrieve a set of field values from a `CplState` variable.
During this exchange, the coupler manages ancillary operations such as 
regridding, unit conversions, or filtering.

## Coupler Object API

```@docs
    ClimaCoupler.CouplerState
    ClimaCoupler.coupler_add_field!
    ClimaCoupler.coupler_add_map!
    ClimaCoupler.coupler_put!
    ClimaCoupler.coupler_get
    ClimaCoupler.coupler_get!
```
