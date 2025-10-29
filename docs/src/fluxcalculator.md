# FluxCalculator

This module contains the infrastructure to compute turbulent fluxes.

## How are fluxes computed?

The key function that computes surface fluxes is
[`FluxCalculator.turbulent_fluxes!`](@ref). This function computes
turbulent fluxes and ancillary quantities, such as the Obukhov length, using
[SurfaceFluxes.jl](https://github.com/CliMA/SurfaceFluxes.jl). Generally, this
function is called at the end of each coupling step.

All the quantities computed in `turbulent_fluxes!` are calculated
separately for each surface model using the
[`FluxCalculator.compute_surface_fluxes!`](@ref) function. This function can be
extended by component models if they need specific type of flux calculation, and
a default is provided for models that can use the standard flux calculation.

The default method of [`FluxCalculator.compute_surface_fluxes!`](@ref), in turn,
calls [`FluxCalculator.get_surface_fluxes`](@ref). This function uses a thermal
state obtained by using the model surface temperature, extrapolates atmospheric
density adiabatically to the surface, and with the surface humidity (if
available, if not, assuming a saturation specific humidity for liquid phase).
`compute_surface_fluxes!` also updates the component internal fluxes fields via
[`FluxCalculator.update_turbulent_fluxes!`](@ref), and adds the area-weighted
contribution from this component model to the `CoupledSimulation` fluxes fields.

Any extension of `FluxCalculator.compute_surface_fluxes!` for a particular
surface model is also expected to update both the component models' internal
fluxes and the CoupledSimulation object's fluxes fields.

[`FluxCalculator.compute_surface_fluxes!`](@ref) sets:
- the flux of momenta, `F_turb_ρτxz`, `F_turb_ρτyz`;
- the flux of energy due to latent heat, `F_lh`;
- the flux of energy due to sensible heat, `F_sh`;
- the flux of moisture, `F_turb_moisture`;
- the Obukhov length, `L_MO`;
- the buoyancy flux, `buoyancy_flux`;
- the roughness lengths for momentum and buoyancy, `z0m` and `z0b`;
- the evaporation scaling factor, `beta`,
- the frictional velocity `ustar`.

!!! note

    [`FluxCalculator.compute_surface_fluxes!`](@ref) always returns the area weighted sum, even if this is not necessarily the most meaningful operation for a given quantity (e.g., for the Obukhov length). This can be improved in the future, if you know how, please open an issue.

Note also that [`FluxCalculator.turbulent_fluxes!`](@ref) only
computes turbulent fluxes, not radiative fluxes. Currently, these are computed
within the atmospheric model.

## FluxCalculator API

```@docs
    ClimaCoupler.FluxCalculator.turbulent_fluxes!
    ClimaCoupler.FluxCalculator.compute_surface_fluxes!
    ClimaCoupler.FluxCalculator.get_surface_fluxes
    ClimaCoupler.FluxCalculator.update_turbulent_fluxes!
```
