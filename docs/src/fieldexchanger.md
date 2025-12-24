# FieldExchanger

This module contains general functions for the exchange of fields between the atmospheric and surface component models.

The main function of the `FieldExchanger` is `exchange!`, which does the following:
- imports atmospheric fields into the coupler fields, via `import_atmos_fields!`
- imports area fraction-weighted surface fields into the coupler fields, via `import_combined_surface_fields!`
- updates all component models with the newly-updated coupler fields, via `update_model_sims!`

The specific fields that are exchanged depend on the requirements of the component models:

The fields imported from the atmosphere to the coupler are specified by extending `FieldExchanger.import_atmos_fields!`
The default `import_atmos_fields!` imports radiative fluxes, liquid precipitation, and snow precipitation.

The fields of a component model that get updated by the coupler are specified by extending `FieldExchanger.update_sim!`.
The default `update_sim!` for an atmosphere model updates the direct and diffuse surface albedos,
the surface temperature, and the turbulent fluxes.
The default `update_sim!` for a surface model updates the air density, radiative fluxes,
liquid precipitation, and snow precipitation.
These updates are done via the `update_field!` function, which must be extended for the
particular variable and component model.
If an `update_field!` function is not defined for a particular component model, it will be ignored.

Note that turbulent fluxes are not updated in `update_sim!`, but rather via
`FluxCalculator.update_turbulent_fluxes!`, where fluxes are computed between
the atmosphere and each surface model.

## FieldExchanger API

```@docs
    FieldExchanger.exchange!
    FieldExchanger.update_sim!
    FieldExchanger.step_model_sims!
    FieldExchanger.update_surface_fractions!
    FieldExchanger.set_caches!
```

## FieldExchanger Internal Functions

```@docs
    FieldExchanger.combine_surfaces!
    FieldExchanger.resolve_area_fractions!
    FieldExchanger.import_atmos_fields!
```
