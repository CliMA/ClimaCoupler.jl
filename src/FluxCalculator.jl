"""
    FluxCalculator

This modules contains abstract types and functions to calculate surface fluxes in the coupler,
or to call flux calculating functions from the component models.
"""
module FluxCalculator

import LazyBroadcast: @lazy
import StaticArrays
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import ClimaCore as CC
import NVTX
import ..Interfacer, ..Utilities

export turbulent_fluxes!,
    get_surface_params,
    reset_fluxes!,
    update_turbulent_fluxes!,
    compute_surface_fluxes!,
    ocean_seaice_fluxes!,
    FluxAccumulator,
    accumulate!,
    push_and_reset!,
    push_ready_accumulators!,
    reset!

"""
    reset_fluxes!(cs::Interfacer.CoupledSimulation)
    reset_fluxes!(sim::Interfacer.AbstractComponentSimulation)

Reset any internally accumulated surface flux tendencies on each component model in `cs`
(or on a single `sim`) at the start of a coupler timestep, before flux contributions are
added in `FieldExchanger.exchange!`/`FieldExchanger.update_sim!`,
[`FluxCalculator.update_turbulent_fluxes!`](@ref), and/or
[`FluxCalculator.ocean_seaice_fluxes!`](@ref).

The `CoupledSimulation` method is intended to be called explicitly from the top-level
coupler loop (`SimCoordinator.step!`) between `update_surface_fractions!` and `exchange!`,
so the reset ordering is visible at the call site rather than hidden inside the exchange.

By default, the per-`sim` method is a no-op. Component models that accumulate fluxes across
multiple coupler updates should extend that method.
"""
function reset_fluxes!(cs::Interfacer.CoupledSimulation)
    for sim in cs.model_sims
        reset_fluxes!(sim)
    end
    return nothing
end

reset_fluxes!(sim::Interfacer.AbstractComponentSimulation) = nothing

function turbulent_fluxes!(cs::Interfacer.CoupledSimulation)
    turbulent_fluxes!(cs.fields, cs.model_sims, cs.thermo_params, cs.flux_accumulators)
    push_ready_accumulators!(cs.model_sims, cs.flux_accumulators, cs.t[] + cs.Δt_cpl)
    return nothing
end

"""
    turbulent_fluxes!(cs::CoupledSimulation)
    turbulent_fluxes!(fields, model_sims, thermo_params, flux_accumulators = (;))

Compute turbulent fluxes and associated quantities. Store the results in `fields` as
area-weighted sums.

This function uses `SurfaceFluxes.jl` under the hood.

For any surface present in `flux_accumulators`, the per-surface flux is added to
that surface's `FluxAccumulator` instead of being pushed directly to the surface
via `update_turbulent_fluxes!`. The area-weighted combined `cs.fields.F_*` fields
(which the atmosphere reads) are updated every call regardless.

Args:
- `csf`: [Field of NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `AbstractComponentSimulation`s.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.
- `flux_accumulators`: [NamedTuple] of `FluxAccumulator`s keyed by surface
    simulation name (the same keys as `model_sims`). Entries are only present
    for slow explicit surfaces; the default empty NamedTuple disables
    accumulation entirely.

(NB: Radiation surface fluxes are calculated by the atmosphere.)
"""
function turbulent_fluxes!(csf, model_sims, thermo_params, flux_accumulators = (;))
    boundary_space = axes(csf)
    atmos_sim = model_sims.atmos_sim
    FT = CC.Spaces.undertype(boundary_space)

    # Reset the coupler fields will compute. We need to do this because we will compute
    # area-weighted averages
    for p in (:F_turb_ρτxz, :F_turb_ρτyz, :F_lh, :F_sh, :F_turb_moisture)
        fill!(getproperty(csf, p), 0)
    end

    # Compute the surface fluxes for each surface model and add them to `csf`.
    for (name, sim) in pairs(model_sims)
        # If the simulation is an implicit flux simulation, the fluxes are computed in the
        # component model's `step!` function, so we don't need to compute them here.
        sim isa Interfacer.AbstractImplicitFluxSimulation && continue
        # `accumulator` is `nothing` for fast surfaces (and for atmos, which is
        # also iterated but whose `compute_surface_fluxes!` is a no-op).
        accumulator = get(flux_accumulators, name, nothing)
        compute_surface_fluxes!(csf, sim, atmos_sim, thermo_params, accumulator)
    end
    return nothing
end

"""
    FluxAccumulator{F}

Time-accumulator for the five turbulent flux fields between a surface model and the atmosphere,
used for slow surfaces (having `dt > Δt_cpl`) whose fluxes are computed explicitly.
The fields hold the running sum of per-surface (not area-weighted) turbulent fluxes
computed at each coupling step, and `n_steps` counts the number of contributions
added since the last reset.

Used by `turbulent_fluxes!` to push the time-averaged flux to a
slow surface just before it steps. Allocated only for slow explicit
surfaces; fast surfaces and implicit-flux surfaces do not have an accumulator.
"""
struct FluxAccumulator{F <: CC.Fields.Field}
    F_lh::F
    F_sh::F
    F_turb_moisture::F
    F_turb_ρτxz::F
    F_turb_ρτyz::F
    n_steps::Base.RefValue{Int}
end

"""
    FluxAccumulator(boundary_space)

Construct a zero-initialized `FluxAccumulator` whose fields live on the coupler
boundary space (no regridding is needed during accumulation).
"""
function FluxAccumulator(boundary_space)
    return FluxAccumulator(
        CC.Fields.zeros(boundary_space),
        CC.Fields.zeros(boundary_space),
        CC.Fields.zeros(boundary_space),
        CC.Fields.zeros(boundary_space),
        CC.Fields.zeros(boundary_space),
        Ref(0),
    )
end

"""
    accumulate!(acc::FluxAccumulator, fields)

Add the per-surface turbulent fluxes in the NamedTuple `fields` (`F_lh`,
`F_sh`, `F_turb_moisture`, `F_turb_ρτxz`, `F_turb_ρτyz`) into the accumulator
and increment `n_steps`. Called once per coupling step from
`update_flux_fields!` for each slow surface, in place of the direct
`update_turbulent_fluxes!` push.
"""
function accumulate!(acc::FluxAccumulator, fields)
    @. acc.F_lh += fields.F_lh
    @. acc.F_sh += fields.F_sh
    @. acc.F_turb_moisture += fields.F_turb_moisture
    @. acc.F_turb_ρτxz += fields.F_turb_ρτxz
    @. acc.F_turb_ρτyz += fields.F_turb_ρτyz
    acc.n_steps[] += 1
    return nothing
end

"""
    push_and_reset!(sim, acc::FluxAccumulator)

Compute the time-averaged flux (dividing the accumulator fields in-place by
`n_steps`), push it to the surface via `update_turbulent_fluxes!(sim, ...)`,
then zero the accumulator. A no-op if `n_steps` is zero.
"""
function push_and_reset!(sim, acc::FluxAccumulator)
    n = acc.n_steps[]
    iszero(n) && return nothing
    # In-place division avoids allocating
    @. acc.F_lh /= n
    @. acc.F_sh /= n
    @. acc.F_turb_moisture /= n
    @. acc.F_turb_ρτxz /= n
    @. acc.F_turb_ρτyz /= n
    fields = (;
        F_lh = acc.F_lh,
        F_sh = acc.F_sh,
        F_turb_moisture = acc.F_turb_moisture,
        F_turb_ρτxz = acc.F_turb_ρτxz,
        F_turb_ρτyz = acc.F_turb_ρτyz,
    )
    update_turbulent_fluxes!(sim, fields)
    reset!(acc)
    return nothing
end

"""
    reset!(acc::FluxAccumulator)

Zero all accumulator fields and reset the step counter. Called by
`push_and_reset!` after pushing the averaged flux to the surface.
"""
function reset!(acc::FluxAccumulator)
    fill!(acc.F_lh, 0)
    fill!(acc.F_sh, 0)
    fill!(acc.F_turb_moisture, 0)
    fill!(acc.F_turb_ρτxz, 0)
    fill!(acc.F_turb_ρτyz, 0)
    acc.n_steps[] = 0
    return nothing
end

"""
    push_ready_accumulators!(model_sims, flux_accumulators, t_next; force = false)

For each surface model present in `flux_accumulators`, check
whether the surface will step at time `t_next`. If so, compute the
time-averaged flux from the accumulator and write it to the surface boundary
conditions via `push_and_reset!`.

Called by `turbulent_fluxes!(cs)` with `t_next = cs.t[] + cs.Δt_cpl`
immediately after accumulation.

Pass `force = true` to skip the `will_step` check and unconditionally push
every non-empty accumulator to its surface. This is used during initialization
to populate slow-surface BCs before the simulation starts.
"""
function push_ready_accumulators!(
    model_sims,
    flux_accumulators,
    t_next;
    force::Bool = false,
)
    for name in keys(flux_accumulators)
        if force || Interfacer.will_step(model_sims[name], t_next)
            push_and_reset!(model_sims[name], flux_accumulators[name])
        end
    end
    return nothing
end

"""
    get_surface_fluxes(surface_fluxes_params, u_int, T_int, ..., config,
                       update_T_sfc_cb, update_q_vap_sfc;
                       roughness_inputs, scheme, solver_opts, flux_specs)

Uses SurfaceFluxes.jl to calculate turbulent surface fluxes.
Fluxes are computed over the entire surface, even where the relevant surface model is not present.

`update_T_sfc_cb` and `update_q_vap_sfc` are positional so they participate in
broadcast (they may vary per element).  The remaining extended-API arguments
(`roughness_inputs`, `scheme`, `solver_opts`, `flux_specs`) are keyword arguments.
"""
function get_surface_fluxes(
    surface_fluxes_params::SF.Parameters.SurfaceFluxesParameters,
    u_int,
    T_int,
    q_tot_int,
    q_liq_int,
    q_ice_int,
    ρ_int,
    h_int,
    u_sfc,
    T_sfc,
    q_vap_sfc,
    h_sfc,
    d,
    config,
    update_T_sfc_cb = nothing,
    update_q_vap_sfc = nothing;
    roughness_inputs = nothing,
    scheme = SF.PointValueScheme(),
    solver_opts = nothing,
    flux_specs = nothing,
)
    Φ_sfc = SFP.grav(surface_fluxes_params) * h_sfc
    Δz = h_int - h_sfc

    outputs = SF.surface_fluxes(
        surface_fluxes_params,
        T_int,
        q_tot_int,
        q_liq_int,
        q_ice_int,
        ρ_int,
        T_sfc,
        q_vap_sfc,
        Φ_sfc,
        Δz,
        d,
        u_int,
        u_sfc,
        roughness_inputs,
        config,
        scheme,
        solver_opts,
        flux_specs,
        update_T_sfc_cb,
        update_q_vap_sfc,
    )

    (; shf, lhf, evaporation, ρτxz, ρτyz, T_sfc, q_vap_sfc) = outputs

    return (;
        F_turb_ρτxz = ρτxz,
        F_turb_ρτyz = ρτyz,
        F_sh = shf,
        F_lh = lhf,
        F_turb_moisture = evaporation,
        T_sfc_new = T_sfc,
    )
end

"""
    get_surface_params(atmos_sim::Interfacer.AbstractAtmosSimulation)

Returns the surface parameters of type `SF.Parameters.SurfaceFluxesParameters`.

TODO: in the future this may not need to depend on the atmos sim, but
here retaining the dependency until we know how EDMF boundary conditions will
be handled (for consistency of parameters).
"""
function get_surface_params(atmos_sim::Interfacer.AbstractAtmosSimulation)
    return error(
        "get_surface_params is required to be dispatched on $(nameof(atmos_sim)), but no method defined",
    )
end

"""
    update_turbulent_fluxes!(sim::Interfacer.AbstractComponentSimulation, fields::NamedTuple)

Updates the fluxes in the simulation `sim` with the fluxes in `fields`.

For surface models, this should be the fluxes computed between the surface model and the atmosphere.
For atmosphere models, this should be the area-weighted sum of fluxes across all surface models.
"""
function update_turbulent_fluxes!(sim::Interfacer.AbstractComponentSimulation, fields)
    return error(
        "update_turbulent_fluxes! is required to be dispatched on $(nameof(sim)), but no method defined",
    )
end

update_turbulent_fluxes!(sim::Interfacer.AbstractSurfaceStub, fields::NamedTuple) = nothing

"""
    compute_surface_fluxes!(csf, sim, atmos_sim, thermo_params)

This function computes surface fluxes between the input component model
simulation and the atmosphere.

Update the input coupler surface fields `csf` in-place with the computed fluxes
for this model. These are then summed using area-weighting across all surface
models to get the total fluxes.

Since the fluxes are computed between the input model and the atmosphere, this
function does nothing if called on an atmosphere model simulation.

The function for AbstractImplicitFluxSimulation is a placeholder that does nothing. Currently,
the only AbstractImplicitFluxSimulation is ClimaLandSimulation, for which compute_surface_fluxes!
is defined in the component model. We can extend this function for other AbstractImplicitFluxSimulation
in the future.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_lh`, `F_sh`, `F_turb_moisture`.
- `sim`: [Interfacer.AbstractComponentSimulation] the surface simulation to compute fluxes for.
- `atmos_sim`: [Interfacer.AbstractAtmosSimulation] the atmosphere simulation to compute fluxes with.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.

The roughness model is obtained from the simulation via `get_field(sim, Val(:roughness_model))`.
Ocean simulations return `:coare3`, while land and ice simulations return `:constant` (the default).
"""
function compute_surface_fluxes!(
    csf,
    sim::Interfacer.AbstractAtmosSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    thermo_params,
    accumulator = nothing,
)
    # do nothing for atmos model
    return nothing
end

function compute_surface_fluxes!(
    csf,
    sim::Interfacer.AbstractImplicitFluxSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    thermo_params,
    accumulator = nothing,
)
    # do nothing for implicit flux surface model
    return nothing
end

NVTX.@annotate function compute_surface_fluxes!(
    csf,
    sim::Interfacer.AbstractSurfaceSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    thermo_params,
    accumulator = nothing,
)
    boundary_space = axes(csf)
    FT = CC.Spaces.undertype(boundary_space)
    surface_fluxes_params = FluxCalculator.get_surface_params(atmos_sim)

    # Atmosphere fields are stored in coupler fields so we only regrid them once per timestep
    # `_int` refers to atmos state of center level 1
    uv_int = StaticArrays.SVector.(csf.u_int, csf.v_int)

    # compute surface humidity from the surface temperature, surface density, and phase
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:surface_temperature))
    T_sfc = csf.scalar_temp1

    # TODO: This is not accurate - we shouldn't assume condensate is 0.
    ρ_sfc = SF.surface_density.(
        surface_fluxes_params,
        csf.T_atmos,
        csf.ρ_atmos,
        T_sfc,
        csf.height_int .- csf.height_sfc,
        csf.q_tot_atmos,
        0, # q_liq
        0, # q_ice
    )

    csf.scalar_temp2 .= TD.q_vap_saturation.(thermo_params, T_sfc, ρ_sfc, 0, 0)
    q_sfc = csf.scalar_temp2

    # Set SurfaceFluxConfig containing models for roughness and gustiness
    roughness_params = get_roughness_params(csf, sim)
    gustiness_spec = SF.ConstantGustinessSpec(FT(1))
    config = @lazy SF.SurfaceFluxConfig.(roughness_params, gustiness_spec)

    # Set surface velocity to zero for now
    uv_sfc = @lazy @. (uv_int * FT(0))
    fluxes = FluxCalculator.get_surface_fluxes.(
        surface_fluxes_params,
        uv_int,
        csf.T_atmos,
        csf.q_tot_atmos,
        csf.q_liq_atmos,
        csf.q_ice_atmos,
        csf.ρ_atmos,
        csf.height_int, # h_int
        uv_sfc,
        T_sfc,
        q_sfc,
        csf.height_sfc, # h_sfc
        FT(0), # d
        config,
    )

    # Update the coupler fields and surface simulation with the fluxes
    update_flux_fields!(csf, sim, fluxes, accumulator)
    return nothing
end

"""
    update_flux_fields!(csf, sim::Interfacer.AbstractSurfaceSimulation, fluxes, accumulator = nothing)

Update the surface simulation `sim` with the values stored in `fluxes`,
without any area-weighting. Then, update the coupler fields `csf` with
the area-weighted fluxes.

If `accumulator` is a `FluxAccumulator` (i.e. `sim` is a slow surface), the
per-surface fluxes are added to the accumulator instead of being pushed
directly to the surface's boundary conditions. The
area-weighted contribution to `csf.F_*` (which the atmosphere reads each
coupling step) is updated unconditionally.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields:
    `F_turb_ρτxz`, `F_turb_ρτyz`, `F_lh`, `F_sh`, `F_turb_moisture`.
- `sim`: [Interfacer.AbstractComponentSimulation] the surface simulation to update.
- `fluxes`: [NamedTuple] containing the fluxes to update the surface simulation with.
- `accumulator`: optional `FluxAccumulator` for slow surfaces; `nothing` (default)
    updates the surface immediately.
"""
function update_flux_fields!(
    csf,
    sim::Interfacer.AbstractSurfaceSimulation,
    fluxes,
    accumulator = nothing,
)
    (; F_turb_ρτxz, F_turb_ρτyz, F_sh, F_lh, F_turb_moisture) = fluxes
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

    # Zero out fluxes where the area fraction is zero
    # Multiplying by `area_fraction` is not sufficient because the fluxes may
    # be NaN where the area fraction is zero.
    @. F_turb_ρτxz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτxz), F_turb_ρτxz)
    @. F_turb_ρτyz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτyz), F_turb_ρτyz)
    @. F_sh = ifelse(area_fraction ≈ 0, zero(F_sh), F_sh)
    @. F_lh = ifelse(area_fraction ≈ 0, zero(F_lh), F_lh)
    @. F_turb_moisture = ifelse(area_fraction ≈ 0, zero(F_turb_moisture), F_turb_moisture)

    # update the fluxes, which are now area fraction-masked, of this surface model
    fields = (; F_turb_ρτxz, F_turb_ρτyz, F_lh, F_sh, F_turb_moisture)

    if isnothing(accumulator)
        FluxCalculator.update_turbulent_fluxes!(sim, fields)
    else
        accumulate!(accumulator, fields)
    end

    # update fluxes in the coupler fields
    # add the flux contributing from this surface to the coupler field
    # note that the fluxes are area-weighted here when we provide them to the atmosphere
    @. csf.F_turb_ρτxz += F_turb_ρτxz * area_fraction
    @. csf.F_turb_ρτyz += F_turb_ρτyz * area_fraction
    @. csf.F_lh += F_lh * area_fraction
    @. csf.F_sh += F_sh * area_fraction
    @. csf.F_turb_moisture += F_turb_moisture * area_fraction

    return nothing
end

"""
    get_roughness_params(csf, sim)

Return the roughness parameters for the simulation, based on the roughness model.
"""
function get_roughness_params(csf, sim)
    roughness_model = Interfacer.get_field(sim, Val(:roughness_model))
    return get_roughness_params(csf, sim, Val(roughness_model))
end

"""
    get_roughness_params(csf, sim, ::Val{:coare3})

Return COARE3 roughness parameters from the simulation.
"""
get_roughness_params(csf, sim, ::Val{:coare3}) =
    Interfacer.get_field(sim, Val(:coare3_roughness_params))

"""
    get_roughness_params(csf, sim, ::Val{:constant})

Load momentum and buoyancy roughness from the simulation into csf scratch fields
and return element-wise `ConstantRoughnessParams`.
"""
function get_roughness_params(csf, sim, ::Val{:constant})
    Interfacer.get_field!(csf.scalar_temp3, sim, Val(:roughness_momentum))
    z0m = csf.scalar_temp3
    Interfacer.get_field!(csf.scalar_temp4, sim, Val(:roughness_buoyancy))
    z0b = csf.scalar_temp4
    return @lazy SF.ConstantRoughnessParams.(z0m, z0b)
end

"""
    get_roughness_params(csf, sim, roughness_model)

Fallback for unknown roughness model; errors.
"""
get_roughness_params(csf, sim, roughness_model) =
    error("Unknown roughness_model: $roughness_model. Must be :coare3 or :constant")


"""
    ocean_seaice_fluxes!(cs::CoupledSimulation)
    ocean_seaice_fluxes!(ocean_sim, ice_sim)

Compute the fluxes between the ocean and sea ice simulations.
This function does nothing by default - it should be extended
for any ocean and sea ice models that support flux calculations.
"""
function ocean_seaice_fluxes!(cs::Interfacer.CoupledSimulation)
    haskey(cs.model_sims, :ocean_sim) &&
        haskey(cs.model_sims, :ice_sim) &&
        ocean_seaice_fluxes!(cs.model_sims.ocean_sim, cs.model_sims.ice_sim)
    return nothing
end
function ocean_seaice_fluxes!(
    ocean_sim::Union{Interfacer.AbstractOceanSimulation, Interfacer.AbstractSurfaceStub},
    ice_sim::Union{Interfacer.AbstractSeaIceSimulation, Interfacer.AbstractSurfaceStub},
)
    return nothing
end

end # module
