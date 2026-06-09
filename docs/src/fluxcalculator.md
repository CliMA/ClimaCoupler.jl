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

!!! note

    [`FluxCalculator.compute_surface_fluxes!`](@ref) always returns the area weighted sum, even if this is not necessarily the most meaningful operation for a given quantity (e.g., for the Obukhov length). This can be improved in the future, if you know how, please open an issue.

Note also that [`FluxCalculator.turbulent_fluxes!`](@ref) only
computes turbulent fluxes, not radiative fluxes. Currently, these are computed
within the atmospheric model.

## Turbulent flux accumulation for slow surfaces

When a surface model's own timestep is larger than the coupling timestep
(e.g. `dt_ocean > Δt_cpl`), pushing a fresh flux to it every coupling step is wasteful:
the surface only consumes its boundary fluxes when it steps, and the intermediate writes
are typically overwritten by the next coupling step exchange. For correctness and
conservation, it is also better to feed the surface a time-average of the per-coupling-step
fluxes rather than a single instantaneous snapshot.

To support this, the coupler allocates a [`FluxCalculator.FluxAccumulator`](@ref)
for each slow explicit surface (see [`Interfacer.sim_dt`](@ref),
[`Interfacer.will_step`](@ref)), and:

- Each coupling step, [`FluxCalculator.turbulent_fluxes!`](@ref) calls
  [`FluxCalculator.accumulate!`](@ref) on the accumulator instead of writing the
  flux directly to the surface via `update_turbulent_fluxes!`. The area-weighted
  combined fields in the CoupledSimulation `fields` NamedTuple, which get sent to the
  atmosphere, are still updated every call.
- [`FluxCalculator.turbulent_fluxes!`](@ref) also calls
  [`FluxCalculator.push_ready_accumulators!`](@ref), which checks if the surface is
  about to step and if so, divides the accumulator by `n_steps`, writes the time-averaged
  flux to the surface boundary conditions, and zeros the accumulator.

Accumulators are not allocated for fast surfaces (`sim_dt ≤ Δt_cpl`),
`AbstractSurfaceStub`s, or `AbstractImplicitFluxSimulation`s (which compute
their own fluxes inside `step!`).

## FluxCalculator API

```@docs
    FluxCalculator.turbulent_fluxes!
    FluxCalculator.compute_surface_fluxes!
    FluxCalculator.get_surface_fluxes
    FluxCalculator.update_turbulent_fluxes!
    FluxCalculator.update_flux_fields!
    FluxCalculator.get_roughness_params
    FluxCalculator.FluxAccumulator
    FluxCalculator.accumulate!
    FluxCalculator.push_and_reset!
    FluxCalculator.push_ready_accumulators!
    FluxCalculator.reset!
```
