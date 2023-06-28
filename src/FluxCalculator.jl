"""
    FluxCalculator

This modules contains abstract types and functions to calculate surface fluxes in the coupler,
or to call flux calculating functions from the component models.
"""

module FluxCalculator

using ClimaCoupler.Interfacer
using ClimaCore.Fields
export PartitionedComponentModelGrid,
    CombinedAtmosGrid,
    compute_combined_turbulent_fluxes!,
    TurbulentFluxPartition,
    compute_atmos_turbulent_fluxes!,
    calculate_surface_air_density

"""
    TurbulentFluxPartition

Abstract type for flags that denote where and how to calculate tubulent fluxes.
"""
abstract type TurbulentFluxPartition end

"""
    PartitionedComponentModelGrid <: TurbulentFluxPartition

A flag indicating that the turbulent fluxes should be partitioned and calculated
over each surface model and then combined. This is calculated on the coupler grid.
"""
struct PartitionedComponentModelGrid <: TurbulentFluxPartition end

"""
    CombinedAtmosGrid <: TurbulentFluxPartition

A flag indicating that the turbulent fluxes (e.g. sensible and latent heat fluxes,
drag and moisture fluxes) are to be  calculated on the Atmos grid, and saved in Atmos cache.
"""
struct CombinedAtmosGrid <: TurbulentFluxPartition end

"""
    compute_combined_turbulent_fluxes!(model_sims, csf, turbulent_fluxes::TurbulentFluxPartition)

Calls the method(s) which calculate turbulent surface fluxes from combined surface states in coupler fields, `csf`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `csf`: [NamedTuple] containing coupler fields.
- `turbulent_fluxes`: [TurbulentFluxPartition] denotes a flag for turbulent flux calculation.

"""
function compute_combined_turbulent_fluxes!(model_sims, csf, turbulent_fluxes::TurbulentFluxPartition)
    if turbulent_fluxes == CombinedAtmosGrid()
        compute_atmos_turbulent_fluxes!(model_sims.atmos_sim, csf)
    else
        nothing # TODO: may want to add CombinedCouplerGrid
    end
end

"""
    compute_atmos_turbulent_fluxes!(sim::Interfacer.ComponentModelSimulation, csf)

A function to calculate turbulent surface fluxes using the combined surface states.
It is required that a method is defined for the given `sim` and that the fluxes are
saved in that sim's cache. `csf` refers to the coupler fields.

# Arguments
- `sim`: [Interfacer.ComponentModelSimulation] object containing the component model simulation.
- `csf`: [NamedTuple] containing coupler fields.

# Example:

```
function compute_atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)
    atmos_sim.cache.flux .= atmos_sim.c .* (csf.T_S .- atmos_sim.temperature)
end
```

"""
compute_atmos_turbulent_fluxes!(sim::Interfacer.ComponentModelSimulation, _) =
    error("calling flux calculation in " * Interfacer.name(sim) * " but no method defined")


"""
    calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_S::Fields.Field)

Extension for this  to to calculate surface density.
"""
function calculate_surface_air_density(atmos_sim::Interfacer.AtmosModelSimulation, T_S::Fields.Field)
    error("this function is required to be dispatched on" * Interfacer.name(atmos_sim) * ", but no method defined")
end

end # module
