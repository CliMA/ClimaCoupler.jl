# FieldExchanger

This module contains general functions for the exchange of fields between the atmospheric and surface component models.

The `FieldExchanger` needs to populate the coupler with
- atmospheric fields, via the `import_atmos_fields!` function
- average surface properties of each coupler gridpoint, via the `import_combined_surface_fields!` function

The component models are updated by broadcasting the coupler fields, via the `update_model_sims!` function. For an update, this function requires that `update_field!` is defined for the particular variable and component model.
If an `update_field!` function is not defined for a particular component model, it will be ignored.

The specific fields that are exchanged depend on the requirements of the component models:

The fields imported from the atmosphere to the coupler are specified by extending `FieldExchanger.import_atmos_fields!`
The default `import_atmos_fields!` imports radiative fluxes, liquid precipitation, and snow precipitation.

The fields of a component model that get updated by the coupler are specified by extending `FieldExchanger.update_sim!`
The default `update_sim!` for an atmosphere model updates the direct
and diffuse surface albedos, the surface temperature,
and the turbulent fluxes.
The default `update_sim!` for a surface model updates the air density, radiative fluxes,
liquid precipitation, and snow precipitation.

## FieldExchanger API

```@docs
    ClimaCoupler.FieldExchanger.import_atmos_fields!
    ClimaCoupler.FieldExchanger.import_combined_surface_fields!
    ClimaCoupler.FieldExchanger.update_model_sims!
    ClimaCoupler.FieldExchanger.update_sim!
    ClimaCoupler.FieldExchanger.reinit_model_sims!
    ClimaCoupler.FieldExchanger.step_model_sims!
    ClimaCoupler.FieldExchanger.update_surface_fractions!
```

## FieldExchanger Internal Functions

```@docs
    ClimaCoupler.FieldExchanger.combine_surfaces!
    ClimaCoupler.FieldExchanger.dummmy_remap!
```
