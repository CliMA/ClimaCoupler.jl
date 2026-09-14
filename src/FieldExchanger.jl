"""
    FieldExchanger

This modules contains general functions for the exchange of fields between the
atmospheric and surface component models.
"""
module FieldExchanger

import ClimaCore as CC
import .Threads

import ..Interfacer, ..FluxCalculator, ..Utilities

"True when the coupler fields carry the overlap cache (i.e. `overlap_slow_surfaces` is on)."
has_slow_cache(csf, val) = slow_cache_name(val) in propertynames(csf)

"Coupler field parking the slow surfaces' contribution to each blended quantity."
slow_cache_name(::Val{:emissivity}) = :slow_emissivity
slow_cache_name(::Val{:surface_temperature}) = :slow_LW_up
slow_cache_name(::Val{:surface_direct_albedo}) = :slow_direct_albedo
slow_cache_name(::Val{:surface_diffuse_albedo}) = :slow_diffuse_albedo

export update_sim!,
    update_model_sims!,
    step_model_sims!,
    exchange!,
    set_caches!,
    step_slow_sims!,
    slow_step_target,
    slow_launch_target,
    slow_progress_snapshot,
    slow_sim_dt,
    slow_window_steps,
    slow_step_boundary,
    launch_slow_sims!,
    wait_slow_sims!,
    slow_step_in_flight,
    update_surface_fractions!,
    resolve_area_fractions!,
    align_surface_fractions!

"""
    update_surface_fractions!(cs::Interfacer.CoupledSimulation)

Updates dynamically changing area fractions.
Maintains the invariant that the sum of area fractions is 1 at all points.
Area fractions are expected to be defined on the boundary space of the coupled simulation,
since they are used by the coupler.

If a surface model is not present, the area fraction is set to 0.

# Arguments
- `cs`: [Interfacer.CoupledSimulation] containing area fraction information.
"""
function update_surface_fractions!(
    cs::Interfacer.CoupledSimulation;
    slow_frozen::Bool = false,
)
    # Area fractions are derived from the sea ice concentration and the ocean
    # wet mask, so they cannot be recomputed while those models are stepping.
    # They are also unchanged over the window, since neither model advances the
    # coupler's view of its state until the step is joined.
    slow_frozen && return nothing

    # An ocean model may provide its own authoritative surface fractions
    # (e.g. derived from its bathymetric wet mask); if it does, skip the
    # land-fraction-based derivation below.
    if haskey(cs.model_sims, :ocean_sim) &&
       align_surface_fractions!(cs.model_sims.ocean_sim, cs)
        return nothing
    end
    _update_surface_fractions_from_land_fraction!(cs)
    return nothing
end

"""
    align_surface_fractions!(ocean_sim, cs::Interfacer.CoupledSimulation) -> Bool

Give the ocean model the opportunity to set all surface area fractions from its
own representation of the land/sea distribution (e.g. a fraction derived from
the ocean bathymetry via an exchange grid).

Return `true` if the fractions were updated (in which case the land-fraction-based update in
[`update_surface_fractions!`](@ref) is skipped), `false` otherwise. The default
implementation returns `false`; ocean models can extend this method.

Implementations must maintain the invariant that the land, ice, and ocean area
fractions sum to 1 at every point of the boundary space.
"""
function align_surface_fractions!(ocean_sim, cs::Interfacer.CoupledSimulation)
    return false
end

function _update_surface_fractions_from_land_fraction!(cs::Interfacer.CoupledSimulation)
    FT = CC.Spaces.undertype(Interfacer.boundary_space(cs))

    # land fraction is static
    if haskey(cs.model_sims, :land_sim)
        land_fraction = Interfacer.get_field(cs.model_sims.land_sim, Val(:area_fraction))
    else
        cs.fields.scalar_temp1 .= 0
        land_fraction = cs.fields.scalar_temp1
    end

    # ice and ocean fractions are dynamic
    if haskey(cs.model_sims, :ice_sim)
        ice_sim = cs.model_sims.ice_sim
        Interfacer.get_field!(cs.fields.scalar_temp2, ice_sim, Val(:ice_concentration))
        ice_concentration = cs.fields.scalar_temp2

        # max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
        Interfacer.update_field!(
            ice_sim,
            Val(:area_fraction),
            max.(min.(ice_concentration, FT(1) .- land_fraction), FT(0)),
        )
        ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))
    else
        cs.fields.scalar_temp2 .= 0
        ice_fraction = cs.fields.scalar_temp2
    end

    if haskey(cs.model_sims, :ocean_sim)
        ocean_sim = cs.model_sims.ocean_sim
        Interfacer.update_field!(
            ocean_sim,
            Val(:area_fraction),
            max.(FT(1) .- ice_fraction .- land_fraction, FT(0)),
        )
        ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction))

        # Apply any additional constraints on the ocean and ice fractions if necessary
        if haskey(cs.model_sims, :ice_sim)
            resolve_area_fractions!(ocean_sim, cs.model_sims.ice_sim, land_fraction)
        end
    else
        cs.fields.scalar_temp3 .= 0
        ocean_fraction = cs.fields.scalar_temp3
    end

    # update the ice and ocean area fraction coupler fields (land is static)
    cs.fields.ice_area_fraction .= ice_fraction
    cs.fields.ocean_area_fraction .= ocean_fraction

    # check that the sum of area fractions is 1
    @assert minimum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
    @assert maximum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
end

"""
    resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)

Ensure that the ocean and ice fractions are consistent with each other.
For most ocean and ice models, this does nothing since the ocean fraction is
defined as `1 - ice_fraction - land_fraction`. However, some models may have
additional constraints on the ice and ocean fractions that need to be enforced.
This function can be extended for such models.
"""
function resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)
    return nothing
end

"""
    import_atmos_fields!(csf, model_sims)

Update the coupler with quantities from the  atmosphere model. By default, this
updates the coupler fields for quantities required for turbulent flux calculations,
radiative fluxes, and precipitation.
This function should be extended for any model that requires additional fields
from the atmosphere.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `AbstractComponentSimulation`s.
"""
function import_atmos_fields!(csf, model_sims)
    # get atmosphere properties used for flux calculations
    Interfacer.get_field!(csf.T_atmos, model_sims.atmos_sim, Val(:air_temperature))
    Interfacer.get_field!(
        csf.q_tot_atmos,
        model_sims.atmos_sim,
        Val(:total_specific_humidity),
    )
    Interfacer.get_field!(
        csf.q_liq_atmos,
        model_sims.atmos_sim,
        Val(:liquid_specific_humidity),
    )
    Interfacer.get_field!(
        csf.q_ice_atmos,
        model_sims.atmos_sim,
        Val(:ice_specific_humidity),
    )
    Interfacer.get_field!(csf.ρ_atmos, model_sims.atmos_sim, Val(:air_density))
    Interfacer.get_field!(csf.u_int, model_sims.atmos_sim, Val(:u_int))
    Interfacer.get_field!(csf.v_int, model_sims.atmos_sim, Val(:v_int))

    # radiative fluxes
    Interfacer.get_field!(csf.SW_d, model_sims.atmos_sim, Val(:SW_d))
    Interfacer.get_field!(csf.LW_d, model_sims.atmos_sim, Val(:LW_d))

    # precipitation
    Interfacer.get_field!(csf.P_liq, model_sims.atmos_sim, Val(:liquid_precipitation))
    Interfacer.get_field!(csf.P_snow, model_sims.atmos_sim, Val(:snow_precipitation))

    for sim in model_sims
        import_atmos_fields!(csf, sim, model_sims.atmos_sim)
    end
end

"""
    import_atmos_fields!(csf, ::Interfacer.AbstractComponentSimulation, atmos_sim)

Updates the coupler simulation fields with atmospheric fluxes from the atmosphere simulation.
This function should be extended for any surface model that requires additional fields
from the atmosphere. Any fields added in a method of this function should also be added
in the corresponding method of `Interfacer.add_coupler_fields!`. The combination of these
two functions defines any extra atmosphere fields provided to the surface.
"""
import_atmos_fields!(csf, ::Interfacer.AbstractComponentSimulation, atmos_sim) = nothing

"""
    import_combined_surface_fields!(csf, model_sims)

Updates the coupler with the surface properties. The `Interfacer.get_field`
functions for (`:emissivity`, `:surface_temperature`, `:surface_direct_albedo`,
`:surface_diffuse_albedo`) need to be specified for each surface model.

Note: Not all surface fields are imported here. Some quantities are retrieved
from each surface model when surface fluxes are computed, in `compute_surface_fluxes!`.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `AbstractComponentSimulation`s.
"""
function import_combined_surface_fields!(csf, model_sims; slow_frozen::Bool = false)
    combine_surfaces!(csf, model_sims, Val(:emissivity); slow_frozen)
    combine_surfaces!(csf, model_sims, Val(:surface_temperature); slow_frozen)
    combine_surfaces!(csf, model_sims, Val(:surface_direct_albedo); slow_frozen)
    combine_surfaces!(csf, model_sims, Val(:surface_diffuse_albedo); slow_frozen)
    return nothing
end

"""
    import_static_fields!(csf, model_sims)

Import static fields into the coupler fields.
This is used to import fields that are not updated during a simulation,
so it is only called at initialization.

Fields imported here are:
- the bottom cell center and face heights of the atmosphere

Any fields imported into the coupler fields here that need to be sent to
component models should be updated in the `set_cache!` function
of the receiving component model.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `AbstractComponentSimulation`s.
"""
function import_static_fields!(csf, model_sims)
    Interfacer.get_field!(csf.height_int, model_sims.atmos_sim, Val(:height_int))
    Interfacer.get_field!(csf.height_sfc, model_sims.atmos_sim, Val(:height_sfc))
    csf.height_delta .= Interfacer.get_atmos_height_delta(csf.height_int, csf.height_sfc)

    return nothing
end

"""
    update_sim!(atmos_sim::Interfacer.AbstractAtmosSimulation, csf)

Updates the atmosphere's fields for surface direct and diffuse albedos, emissivity,and temperature.

# Arguments
- `atmos_sim`: [Interfacer.AbstractAtmosSimulation] containing an atmospheric model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(atmos_sim::Interfacer.AbstractAtmosSimulation, csf)
    Interfacer.update_field!(
        atmos_sim,
        Val(:surface_direct_albedo),
        csf.surface_direct_albedo,
    )
    Interfacer.update_field!(
        atmos_sim,
        Val(:surface_diffuse_albedo),
        csf.surface_diffuse_albedo,
    )
    Interfacer.update_field!(atmos_sim, Val(:emissivity), csf.emissivity)
    Interfacer.update_field!(atmos_sim, Val(:surface_temperature), csf.T_sfc)
    Interfacer.update_field!(atmos_sim, Val(:surface_humidity), csf)
    return nothing
end

"""
    update_sim!(sim::AbstractSurfaceSimulation, csf)

Updates the surface component model cache with the current coupler fields
*besides turbulent fluxes*, which are updated in `update_turbulent_fluxes`.

Note that upwelling longwave and shortwave radiation are not computed here,
and are expected to be computed internally by the surface model.
Some component models extend this function and compute the upwelling longwave
and shortwave radiation in their methods of `update_sim!`.

# Arguments
- `sim`: [Interfacer.AbstractSurfaceSimulation] containing a surface model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(sim::Interfacer.AbstractSurfaceSimulation, csf)
    # radiative fluxes
    Interfacer.update_field!(sim, Val(:SW_d), csf.SW_d)
    Interfacer.update_field!(sim, Val(:LW_d), csf.LW_d)

    # precipitation
    Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)
    return nothing
end

"""
    update_model_sims!(model_sims, csf)

Iterates `update_sim!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `AbstractComponentSimulation`s.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_model_sims!(model_sims, csf; slow_frozen::Bool = false)
    for sim in model_sims
        # `update_sim!` for the ocean zeroes and rebuilds its surface flux
        # fields, which the ocean is integrating with while a slow step is in
        # flight. Leave those sims alone until the step is joined.
        slow_frozen && Interfacer.is_overlapped(sim) && continue
        update_sim!(sim, csf)
    end
end

function Interfacer.step!(
    land_sim::Interfacer.AbstractLandSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    t,
    coupler_fields,
    thermo_params,
)
    # Step the land simulation first
    Interfacer.step!(land_sim, t)

    # Update the atmosphere with the fluxes across all surface models. An explicit-flux
    # land model does not precompute fluxes the way an `AbstractImplicitFluxSimulation`
    # does, but the atmosphere still must not be stepped before this update -- the serial
    # branch of `step_model_sims!` always performs it.
    FluxCalculator.update_turbulent_fluxes!(atmos_sim, coupler_fields)

    # Step the atmosphere model
    Interfacer.step!(atmos_sim, t)
end

function Interfacer.step!(
    implicit_flux_sim::Interfacer.AbstractImplicitFluxSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    t,
    coupler_fields,
    thermo_params,
)
    # Step the implicit flux simulation first
    Interfacer.step!(implicit_flux_sim, t)

    # For an implicit flux simulation, `compute_surface_fluxes!` reads in the precomputed
    # fluxes, and puts them into the coupler fields.
    FluxCalculator.compute_surface_fluxes!(
        coupler_fields,
        implicit_flux_sim,
        atmos_sim,
        thermo_params,
    )

    # Update the atmosphere with the fluxes across all surface models
    # The surface models have already been updated with the fluxes in `compute_surface_fluxes!`,
    # or internally within the step in the case of the integrated land model.
    FluxCalculator.update_turbulent_fluxes!(atmos_sim, coupler_fields)

    # Step the atmosphere model
    Interfacer.step!(atmos_sim, t)
end

"""
    step_model_sims!(model_sims, t, coupler_fields, thermo_params)
    step_model_sims!(cs::CoupledSimulation)

Iterates `step!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `AbstractComponentSimulation`s.
- `t`: [AbstractFloat or ITime] denoting the simulation time.
- `coupler_fields`: [Field of NamedTuple] containing the coupler exchange fields.
- `thermo_params`: thermodynamic parameters.
"""
function step_model_sims!(
    model_sims,
    t,
    coupler_fields,
    thermo_params,
    step_concurrently;
    skip_slow::Bool = false,
)
    if step_concurrently && (haskey(model_sims, :ocean_sim) || haskey(model_sims, :ice_sim))
        # Group 1: land and atmosphere. These are implicitly coupled, so they
        # step sequentially inside a single task.
        land_atmos_group = function ()
            if haskey(model_sims, :land_sim)
                Interfacer.step!(
                    model_sims.land_sim,
                    model_sims.atmos_sim,
                    t,
                    coupler_fields,
                    thermo_params,
                )
            else
                # Same ordering requirement as the grouped land/atmos methods above:
                # the atmosphere needs its turbulent fluxes before it steps.
                FluxCalculator.update_turbulent_fluxes!(model_sims.atmos_sim, coupler_fields)
                Interfacer.step!(model_sims.atmos_sim, t)
            end
        end

        # Group 2: sea ice and ocean. These must NOT be separate tasks. The sea
        # ice model is built holding views into the ocean's live surface fields
        # -- `ocean_surface_velocities` and `ocean_surface_salinity` return
        # `view`s of `ocean.model.velocities.u/v` and `ocean.model.tracers.S`,
        # which are captured in the ice model's `SemiImplicitStress` and
        # `IceWaterThermalEquilibrium`. Stepping them in parallel lets the ice
        # read ocean velocity and salinity while the ocean is writing them.
        # Step the ice first, matching the ordering of the sequential branch
        # below, so the ice sees the ocean state from before the ocean steps.
        ice_ocean_group = () -> step_slow_sims!(model_sims, t)

        if skip_slow
            # An asynchronous ice/ocean step is already in flight (or is being
            # launched separately); advance only the fast group here.
            land_atmos_group()
            return nothing
        end

        if get(ENV, "COUPLER_SEQUENTIAL_GROUPS", "0") in ("1", "true", "TRUE", "yes")
            # Diagnostic mode: identical grouping and ordering to the concurrent
            # path, but with no parallelism between the two groups. Used to
            # isolate whether a discrepancy comes from cross-group concurrency or
            # from the regrouping itself.
            land_atmos_group()
            ice_ocean_group()
        else
            @sync begin
                Threads.@spawn land_atmos_group()
                Threads.@spawn ice_ocean_group()
            end
        end
    else
        # Step all surface models (the ordering doesn't matter here)
        for sim in model_sims
            sim isa Interfacer.AbstractSurfaceSimulation && Interfacer.step!(sim, t)

            # For an implicit flux simulation, `compute_surface_fluxes!` reads in the precomputed
            # fluxes, and puts them into the coupler fields.
            sim isa Interfacer.AbstractImplicitFluxSimulation &&
                FluxCalculator.compute_surface_fluxes!(
                    coupler_fields,
                    sim,
                    model_sims.atmos_sim,
                    thermo_params,
                )
        end

        # Update the atmosphere with the fluxes across all surface models
        # The surface models have already been updated with the fluxes in `compute_surface_fluxes!`,
        # or internally within the step in the case of the integrated land model.
        FluxCalculator.update_turbulent_fluxes!(model_sims.atmos_sim, coupler_fields)

        # Step the atmosphere model
        Interfacer.step!(model_sims.atmos_sim, t)
    end
    return nothing
end

function step_model_sims!(cs::Interfacer.CoupledSimulation; skip_slow::Bool = false)
    step_model_sims!(
        cs.model_sims,
        cs.t[],
        cs.fields,
        cs.thermo_params,
        cs.step_concurrently;
        skip_slow,
    )
end

"""
    step_slow_sims!(model_sims, t)

Advance the overlapped group (sea ice, then ocean) to time `t`.

This is the body handed to the asynchronous task when `overlap_slow_surfaces` is
set. The ordering matches the group in `step_model_sims!`: the sea ice holds
views into the ocean's surface velocity and salinity, so it must step before the
ocean rather than beside it.
"""
function step_slow_sims!(model_sims, t; gather_progress::Bool = false)
    for sim in model_sims
        sim isa Interfacer.AbstractSeaIceSimulation && Interfacer.step!(sim, t)
    end
    for sim in model_sims
        sim isa Interfacer.AbstractOceanSimulation && Interfacer.step!(sim, t)
    end
    # Gather the progress scalars here, where this task still owns the state.
    # Reducing over these fields from the coupling loop would race with the
    # writes above; the reports fall due mid-window, so they read this instead.
    # Only the overlapped path needs them; the per-step concurrent path reports
    # from live state after joining, so it skips these extra reductions.
    gather_progress || return nothing
    names = Symbol[]
    snapshots = Any[]
    for (name, sim) in pairs(model_sims)
        Interfacer.is_overlapped(sim) || continue
        push!(names, name)
        push!(snapshots, Interfacer.progress_snapshot(sim))
    end
    return NamedTuple{Tuple(names)}(Tuple(snapshots))
end

"""
    slow_launch_target(cs)

The model time an about-to-be-launched slow step should advance to: exactly one
slow step, never more.

The caller has already established that a step is due by the next coupling time.
If one is *also* due at the current coupling time, the slow sims are a step
behind the coupler and the target is the current time; otherwise it is the next.
Using the next coupling time unconditionally would ask for two steps at once
whenever the slow timestep equals the coupling timestep.
"""
function slow_launch_target(cs::Interfacer.CoupledSimulation)
    t_now = cs.t[]
    if cs.prime_slow_surfaces
        # Primed, the slow group is already level with the coupler at a window
        # boundary and is about to run one window ahead, so the target is simply
        # one slow step on. `will_step` cannot be used here: it compares against
        # the component clock, which priming has moved.
        return t_now + slow_window_steps(cs) * cs.Δt_cpl
    end
    for sim in cs.model_sims
        Interfacer.is_overlapped(sim) || continue
        # Pass the coupler time through unconverted. Under `use_itime` the
        # component clocks hold `DateTime`s and `will_step` dispatches on
        # `ITime`; coercing to Float64 here selects the Float64 method, which
        # then computes `Float64 - DateTime`.
        Interfacer.will_step(sim, t_now) && return t_now
    end
    return t_now + cs.Δt_cpl
end

"""
    launch_slow_sims!(cs)
    wait_slow_sims!(cs)

Start, and later join, the asynchronous ice/ocean step used by
`overlap_slow_surfaces`. Between the two calls the coupler must not read or
write ice or ocean state; the `slow_frozen` paths through `exchange!`,
`update_surface_fractions!`, `turbulent_fluxes!` and `ocean_seaice_fluxes!`
enforce that, and `step_model_sims!` leaves the slow group alone entirely.
"""
function launch_slow_sims!(cs::Interfacer.CoupledSimulation)
    @assert isnothing(cs.slow_task[]) "a slow-surface step is already in flight"
    model_sims = cs.model_sims
    target = slow_launch_target(cs)
    # The target is recorded alongside the task because the ocean's own clock is
    # being written by that task; scheduling decisions must not read it. The task
    # gets the native time (ITime or Float64) so `step!` dispatches correctly;
    # the recorded copy is in seconds, purely for the scheduling comparison.
    cs.slow_task[] = (;
        task = Threads.@spawn(step_slow_sims!(model_sims, target; gather_progress = true)),
        target = Float64(float(target)),
    )
    return nothing
end

function wait_slow_sims!(cs::Interfacer.CoupledSimulation)
    inflight = cs.slow_task[]
    isnothing(inflight) && return nothing
    # Keep the snapshots the task gathered, so reports falling due during the
    # next window have settled values to print.
    cs.slow_progress[] = fetch(inflight.task)
    cs.slow_task[] = nothing
    return nothing
end

"""
    slow_progress_snapshot(cs, sim_name)

Progress scalars gathered at the end of the most recent completed slow step, or
`nothing` if none has completed yet.
"""
function slow_progress_snapshot(cs::Interfacer.CoupledSimulation, sim_name::Symbol)
    snapshots = cs.slow_progress[]
    (isnothing(snapshots) || !haskey(snapshots, sim_name)) && return nothing
    return snapshots[sim_name]
end

"Model time the in-flight slow step is advancing to, or `nothing`."
slow_step_target(cs::Interfacer.CoupledSimulation) =
    isnothing(cs.slow_task[]) ? nothing : cs.slow_task[].target

"""
    slow_step_boundary(cs)

The coupler time at which the overlapped group's next step should be launched,
i.e. where its own clock currently sits.

Found by asking `will_step` when the group would next step and subtracting one
slow step, rather than reading a component clock directly. `will_step` already
knows how to compare coupler time against each component's clock, which under
`use_itime` is a `DateTime` with a different origin than the coupler's seconds
counter -- a conversion that is easy to get wrong.

Only valid when no slow step is in flight, so it is called at construction.
"""
function slow_step_boundary(cs::Interfacer.CoupledSimulation)
    k = slow_window_steps(cs)
    t = cs.t[]
    for _ in 0:(2k + 2)
        stepping = any(
            sim -> Interfacer.is_overlapped(sim) && Interfacer.will_step(sim, t),
            values(cs.model_sims),
        )
        stepping && return t - k * cs.Δt_cpl
        t = t + cs.Δt_cpl
    end
    # No overlapped sims, or none that will step: never fire.
    return nothing
end

"""
    slow_window_steps(cs)

Number of coupling steps spanned by one slow step. Reads only immutable model
fields, so it is safe to call while a slow step is in flight.
"""
slow_window_steps(cs::Interfacer.CoupledSimulation) =
    round(Int, slow_sim_dt(cs) / Float64(float(cs.Δt_cpl)))

"Shortest timestep among the overlapped sims. Reads only immutable fields."
function slow_sim_dt(cs::Interfacer.CoupledSimulation)
    dt = Inf
    for sim in cs.model_sims
        Interfacer.is_overlapped(sim) || continue
        dt = min(dt, Interfacer.sim_dt(sim))
    end
    return dt
end

slow_step_in_flight(cs::Interfacer.CoupledSimulation) = !isnothing(cs.slow_task[])

"""
    combine_surfaces!(csf, sims, field_name_val::Val{field_name}) where {field_name}

Sums the surface fields specified by `field_name_val`, weighted by the respective area fractions
of all surface simulations. The result is saved in the coupler field specified by `field_name_val`.

For surface temperature, upward longwave radiation is computed from the temperatures
of each surface, weighted by their area fractions, and then the combined temperature
is computed from the combined upward longwave radiation.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
    Note: For the surface temperature, all coupler fields are passed in a NamedTuple.
- `sims`: [NamedTuple] containing simulations.
- `field_name_val`: [Val] containing the name Symbol of the field to be extracted by the `Interfacer.get_field` functions.

# Example
- `combine_surfaces!(temp_field, cs.model_sims, Val(:emissivity))`
"""
function combine_surfaces!(
    csf,
    sims,
    field_name_val::Val{field_name};
    slow_frozen::Bool = false,
) where {field_name}
    # Extract the coupler field we are updating
    combined_field = getproperty(csf, field_name)
    FT = eltype(combined_field)
    combined_field .= zero(FT)

    # The slow surfaces' contribution is summed separately so it can be carried
    # across an overlapped step, during which their state must not be read.
    slow_sum = nothing
    if has_slow_cache(csf, field_name_val)
        slow_sum = getproperty(csf, slow_cache_name(field_name_val))
        slow_frozen || (slow_sum .= zero(FT))
    end

    for sim in sims
        if sim isa Interfacer.AbstractSurfaceSimulation
            # While a slow step is in flight, skip those sims and use the parked sum
            slow_frozen && Interfacer.is_overlapped(sim) && continue

            # Store the area fraction of this simulation in `scalar_temp` and rename for clarity
            Interfacer.get_field!(csf.scalar_temp1, sim, Val(:area_fraction))
            area_fraction = csf.scalar_temp1

            # Remap the surface field onto a coupler temporary field to avoid allocation
            Interfacer.get_field!(csf.scalar_temp2, sim, field_name_val)
            surface_field = csf.scalar_temp2

            # Zero out the contribution from this surface if the area fraction is zero.
            # Note that multiplying by `area_fraction` is not sufficient in the case of NaNs
            contribution =
                area_fraction .* ifelse.(area_fraction .≈ 0, zero(FT), surface_field)
            if !isnothing(slow_sum) && Interfacer.is_overlapped(sim)
                slow_sum .+= contribution
            else
                combined_field .+= contribution
            end
        end
    end
    isnothing(slow_sum) || (combined_field .+= slow_sum)
    return nothing
end
function combine_surfaces!(
    csf,
    sims,
    val::Val{:surface_temperature};
    slow_frozen::Bool = false,
)
    # extract the coupler fields we need to get the surface temperature
    T_sfc = csf.T_sfc
    emissivity_sfc = csf.emissivity

    FT = eltype(T_sfc)
    T_sfc .= zero(FT)

    # As in the generic method, the slow surfaces' share of the upward longwave
    # sum is kept separately so it survives an overlapped step.
    slow_sum = nothing
    if has_slow_cache(csf, val)
        slow_sum = getproperty(csf, slow_cache_name(val))
        slow_frozen || (slow_sum .= zero(FT))
    end

    for sim in sims
        if sim isa Interfacer.AbstractSurfaceSimulation
            slow_frozen && Interfacer.is_overlapped(sim) && continue
            # Store the area fraction and emissivity of this simulation in temp fields
            Interfacer.get_field!(csf.scalar_temp1, sim, Val(:area_fraction))
            area_fraction = csf.scalar_temp1
            Interfacer.get_field!(csf.scalar_temp2, sim, Val(:emissivity))
            emissivity_sim = csf.scalar_temp2

            # Remap the surface field onto a coupler temporary field to avoid allocation
            Interfacer.get_field!(csf.scalar_temp3, sim, Val(:surface_temperature))
            T_sfc_sim = csf.scalar_temp3

            # Zero out the contribution from this surface if the area fraction is zero.
            # Note that multiplying by `area_fraction` is not sufficient in the case of NaNs
            # Compute upward longwave radiation from surface temperature for this simulation
            contribution =
                area_fraction .*
                ifelse.(area_fraction .≈ 0, zero(FT), emissivity_sim .* T_sfc_sim .^ FT(4))
            if !isnothing(slow_sum) && Interfacer.is_overlapped(sim)
                slow_sum .+= contribution
            else
                T_sfc .+= contribution
            end
        end
    end
    isnothing(slow_sum) || (T_sfc .+= slow_sum)
    # Convert the combined upward longwave radiation into a surface temperature
    @. T_sfc = (T_sfc / emissivity_sfc)^FT(1 / 4)
    return nothing
end

"""
    exchange!(cs::Interfacer.CoupledSimulation)

Exchange fields between the surface and atmosphere models.
This is done in 2 steps:
1. Import the atmosphere fields and surface fields into the coupler.
2. Update the component model simulations with the coupler fields.

The order of these steps is important, as importing the surface fields requires
the atmosphere fields to be updated so that surface humidity can be computed.
"""
function exchange!(cs::Interfacer.CoupledSimulation; slow_frozen::Bool = false)
    # Import the atmosphere fields and surface fields into the coupler
    import_atmos_fields!(cs.fields, cs.model_sims)
    import_combined_surface_fields!(cs.fields, cs.model_sims; slow_frozen)

    # Update the component model simulations with the coupler fields
    update_model_sims!(cs.model_sims, cs.fields; slow_frozen)
    return nothing
end

"""
    set_caches!(cs::Interfacer.CoupledSimulation)

Perform any initialization of the component model caches that cannot be
done before the initial exchange. This is useful in handling cache interdependencies
between component models.

For example, the radiation callback in the atmosphere model needs to be
initialized with the surface temperatures, which are only available after the
initial exchange. The integrated land, in turn, requires its drivers in the
cache to be filled with the initial radiation fluxes, so that it can propagate
these to the rest of its cache (e.g. in canopy radative transfer).

This function can also be used to set exchanged fields that are static over the
simulation, since it is only called at initialization.
"""
function set_caches!(cs::Interfacer.CoupledSimulation)
    Interfacer.set_cache!(cs.model_sims.atmos_sim, cs.fields)
    exchange!(cs)
    for sim in cs.model_sims
        sim isa Interfacer.AbstractSurfaceSimulation &&
            Interfacer.set_cache!(sim, cs.fields)
    end
    return nothing
end

end # module
