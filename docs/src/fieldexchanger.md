# FieldExchanger

This module contains general functions for the exchange of fields between the atmospheric and surface component models.

The `FieldExchanger` needs to populate the coupler with
- atmospheric fields (mostly fluxes), via the `import_atmos_fields!` function
- average surface properties of each coupler gridpoint, via the `import_combined_surface_fields!` function

The component models are updated by broadcasting the coupler fields, via the `update_model_sims!` function. For an update, this function requires that `update_field!` is defined for the particular variable and component model. Currently, we support the:
- `AtmosModelSimulation`: `albedo`, `surface_temperature`
    - if calculating fluxes in the atmospheric model: `roughness_momentum`, `roughness_buoyancy`, `beta`
- `SurfaceModelSimulation`: `air_density`, `turbulent_energy_flux`, `turbulent_moisture_flux`, `radiative_energy_flux_sfc`, `liquid_precipitation`, `snow_precipitation`

If an `update_field!` function is not defined for a particular component model, it will be ignored.

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
