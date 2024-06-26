# FluxCalculator

This modules contains abstract types and functions to calculate surface fluxes in the coupler, or to call flux calculating functions from the component models.

Fluxes over a heterogeneous surface (e.g., from a gridpoint where atmospheric cell is overlying both land and ocean) can be handled in two different ways:
1. **Combined fluxes** (called with `CombinedStateFluxesMOST()`)
  - these are calculated by averaging the surface properties for each gridpoint (e.g., land and ocean temperatures, albedos and roughness lengths are averaged, based on their respective area fractions), so the flux is calculated only once per gridpoint of the grid where we calculate fluxes. This is computationally faster, but it makes the fluxes on surface boundaries more diffuse. Currently, we use this method for calculating radiative fluxes in the atmosphere, and turbulent fluxes in the coupler (on the atmospheric grid). The fluxes are calculated in the atmospheric (host) model's cache, which can be setup to avoid allocating coupler fields.
2. **Partitioned fluxes** (called with `PartitionedStateFluxes()`)
  - these are calculated separately for each surface type. It is then the fluxes (rather than the surface states) that are combined and passed to the atmospheric model as a boundary condition. This method ensures that each surface model receives fluxes that correspond to its state properties, resulting in a more accurate model evolution. However, it is more computationally expensive, and requires more communication between the component models.

## FluxCalculator API

```@docs
    ClimaCoupler.FluxCalculator.TurbulentFluxPartition
    ClimaCoupler.FluxCalculator.PartitionedStateFluxes
    ClimaCoupler.FluxCalculator.CombinedStateFluxesMOST
    ClimaCoupler.FluxCalculator.combined_turbulent_fluxes!
    ClimaCoupler.FluxCalculator.atmos_turbulent_fluxes_bulk!
    ClimaCoupler.FluxCalculator.atmos_turbulent_fluxes_most!
    ClimaCoupler.FluxCalculator.surface_inputs
    ClimaCoupler.FluxCalculator.calculate_surface_air_density
    ClimaCoupler.FluxCalculator.get_surface_params
    ClimaCoupler.FluxCalculator.partitioned_turbulent_fluxes!
    ClimaCoupler.FluxCalculator.differentiate_turbulent_fluxes!
    ClimaCoupler.FluxCalculator.get_surface_fluxes!
    ClimaCoupler.FluxCalculator.update_turbulent_fluxes!
    ClimaCoupler.FluxCalculator.extrapolate_ρ_to_sfc
    ClimaCoupler.FluxCalculator.surface_thermo_state
    ClimaCoupler.FluxCalculator.water_albedo_from_atmosphere!
```

## FieldExchanger Internal Functions

```@docs
    ClimaCoupler.FluxCalculator.get_scheme_properties
```
