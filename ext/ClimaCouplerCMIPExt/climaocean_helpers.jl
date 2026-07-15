"""
    surface_flux(f::OC.AbstractField)

Extract the top boundary conditions for the given field.
"""
function surface_flux(f::OC.AbstractField)
    top_bc = f.boundary_conditions.top
    if top_bc isa OC.BoundaryCondition{<:OC.BoundaryConditions.Flux}
        return top_bc.condition
    else
        return nothing
    end
end

"""
    set_from_extrinsic_vector!(vector, grid, u_cc, v_cc)

Given the extrinsic vector components `u_cc` and `v_cc` as `Center, Center`
fields, rotate them onto the target grid and remap to `Face, Center` and
`Center, Face` fields, respectively.
"""
function set_from_extrinsic_vector!(vector, grid, u_cc, v_cc)
    arch = OC.Architectures.architecture(grid)

    # Rotate vector components onto the grid
    OC.Utils.launch!(arch, grid, :xy, _rotate_vector!, u_cc, v_cc, grid)

    # Fill halo regions with the rotated vector components so we can use them to interpolate
    OC.fill_halo_regions!(u_cc)
    OC.fill_halo_regions!(v_cc)

    # Interpolate the vector components to face/center and center/face respectively
    OC.Utils.launch!(
        arch,
        grid,
        :xy,
        _interpolate_vector!,
        vector.u,
        vector.v,
        grid,
        u_cc,
        v_cc,
    )
    return nothing
end

"""
    _rotate_vector!(τx, τy, grid)

Rotate the velocities from the extrinsic coordinate system to the intrinsic
coordinate system.
"""
@kernel function _rotate_vector!(τx, τy, grid)
    # Use `k = 1` to index into the reduced Fields
    i, j = @index(Global, NTuple)
    # Rotate u, v from extrinsic to intrinsic coordinate system
    τxr, τyr = OC.Operators.intrinsic_vector(i, j, 1, grid, τx, τy)
    @inbounds begin
        τx[i, j, 1] = τxr
        τy[i, j, 1] = τyr
    end
end

"""
    _interpolate_vector!(τx, τy, grid, τx_cc, τy_cc)

Interpolate the input fluxes `τx_cc` and `τy_cc`, which are Center/Center
Fields to Face/Center and Center/Face coordinates, respectively.
"""
@kernel function _interpolate_vector!(τx, τy, grid, τx_cc, τy_cc)
    # Use `k = 1` to index into the reduced Fields
    i, j = @index(Global, NTuple)
    @inbounds begin
        τx[i, j, 1] = OC.Operators.ℑxᶠᵃᵃ(i, j, 1, grid, τx_cc)
        τy[i, j, 1] = OC.Operators.ℑyᵃᶠᵃ(i, j, 1, grid, τy_cc)
    end
end

# Momentum remapping helpers (CT types / basis fields are defined in oceananigans.jl)

"""
    contravariant_to_cartesian!(remapping, ρτxz, ρτyz)

Convert the covariant vector components `ρτxz` and `ρτyz` from the
contravariant basis (as they are output by the surface flux calculation)
to the Cartesian basis, writing into `remapping.temp_uv_vec`. These are now
in an extrinsic coordinate system that can be rotated onto the ocean/sea ice
grid by `_rotate_vector!`.

Uses precomputed `ct12_xz` / `ct12_yz` stored in `remapping` to avoid
allocating basis fields every coupling step.
"""
function contravariant_to_cartesian!(remapping, ρτxz, ρτyz)
    (; temp_uv_vec, ct12_xz, ct12_yz) = remapping
    local_geometry = CC.Fields.local_geometry_field(ρτxz)
    @. temp_uv_vec =
        CC.Geometry.UVVector(ρτxz * ct12_xz + ρτyz * ct12_yz, local_geometry)
    return nothing
end

# Non-allocating ClimaCore -> Oceananigans remap.
function Interfacer.remap!(target_field::OC.Field, source_field::CC.Fields.Field, remapping)
    # Get the index of the top level (surface); Nz=1 for 2D fields
    Nz = size(target_field, 3)
    dst = vec(OC.interior(target_field, :, :, Nz))

    CR.regrid!(dst, remapping.remapper_cc_to_oc, source_field)
    return nothing
end

# Allocating ClimaCore -> Oceananigans remap
function Interfacer.remap(
    target_space::Union{
        OC.OrthogonalSphericalShellGrid,
        OC.ImmersedBoundaryGrid,
        OC.LatitudeLongitudeGrid,
    },
    source_field::CC.Fields.Field,
    remapping,
)
    target_field = OC.Field{OC.Center, OC.Center, Nothing}(target_space)
    Interfacer.remap!(target_field, source_field, remapping)
    return target_field
end

# Non-allocating Oceananigans Field -> ClimaCore remap.
function Interfacer.remap!(target_field::CC.Fields.Field, source_field::OC.Field, remapping)
    # Get the index of the top level (surface); Nz=1 for 2D fields
    Nz = size(source_field, 3)
    src = vec(OC.interior(source_field, :, :, Nz))

    CR.regrid!(target_field, remapping.remapper_oc_to_cc, src)
    return nothing
end

# Allocating Oceananigans Field -> ClimaCore remap
function Interfacer.remap(
    target_space::CC.Spaces.AbstractSpace,
    source_field::OC.Field,
    remapping,
)
    target_field = CC.Fields.zeros(target_space)
    Interfacer.remap!(target_field, source_field, remapping)
    return target_field
end

# Oceananigans operation -> ClimaCore remap.
# Prefer specialized `get_field!` methods (e.g. `:surface_temperature`) that
# evaluate Celsius→Kelvin into remapping scratch. This generic path still
# allocates a temporary `Field` for arbitrary operations.
function Interfacer.remap!(
    target_field::CC.Fields.Field,
    operation::OC.AbstractOperations.AbstractOperation,
    remapping,
)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    Interfacer.remap!(target_field, evaluated_field, remapping)
    return nothing
end

# Allocating Oceananigans operation -> ClimaCore remap
function Interfacer.remap(
    target_space::CC.Spaces.AbstractSpace,
    operation::OC.AbstractOperations.AbstractOperation,
    remapping,
)
    target_field = CC.Fields.zeros(target_space)
    Interfacer.remap!(target_field, operation, remapping)
    return target_field
end

# Extend Interfacer.get_field to allow automatic remapping to the target space
function Interfacer.get_field!(
    target_field,
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    quantity,
)
    Interfacer.remap!(target_field, Interfacer.get_field(sim, quantity), sim.remapping)
    return nothing
end

"""
    Interfacer.get_field!(target_field, sim::OceananigansSimulation, ::Val{:surface_temperature})

Write ocean surface temperature in Kelvin into `target_field` without allocating
a temporary 3D `Field` for the Celsius→Kelvin conversion (see #1377 / #1384).
"""
function Interfacer.get_field!(
    target_field,
    sim::OceananigansSimulation,
    ::Val{:surface_temperature},
)
    T = sim.ocean.model.tracers.T
    Nz = size(T, 3)
    scratch = sim.remapping.scratch_field_oc1
    C_to_K = sim.ocean_properties.C_to_K
    OC.interior(scratch, :, :, 1) .= OC.interior(T, :, :, Nz) .+ C_to_K
    Interfacer.remap!(target_field, scratch, sim.remapping)
    return nothing
end

"""
    Interfacer.get_field!(target_field, sim::ClimaSeaIceSimulation, ::Val{:surface_temperature})

Write ice surface temperature in Kelvin into `target_field` using remapping
scratch so the Celsius→Kelvin broadcast does not allocate each coupling step.
"""
function Interfacer.get_field!(
    target_field,
    sim::ClimaSeaIceSimulation,
    ::Val{:surface_temperature},
)
    T_c = top_thermodynamics(sim).top_surface_temperature
    scratch = sim.remapping.scratch_field_oc1
    C_to_K = sim.ice_properties.C_to_K
    OC.interior(scratch, :, :, 1) .= OC.interior(T_c, :, :, 1) .+ C_to_K
    Interfacer.remap!(target_field, scratch, sim.remapping)
    return nothing
end

"""
    Interfacer.get_field!(target_field, sim::ClimaSeaIceSimulation, ::Val{:internal_temperature})

Write ocean–ice interface temperature in Kelvin into `target_field` using
remapping scratch to avoid per-call operation-field allocation.
"""
function Interfacer.get_field!(
    target_field,
    sim::ClimaSeaIceSimulation,
    ::Val{:internal_temperature},
)
    T_c = sim.ocean_ice_interface.temperature
    scratch = sim.remapping.scratch_field_oc1
    C_to_K = sim.ice_properties.C_to_K
    OC.interior(scratch, :, :, 1) .= OC.interior(T_c, :, :, 1) .+ C_to_K
    Interfacer.remap!(target_field, scratch, sim.remapping)
    return nothing
end

# TODO see if we can remove this allocating version
function Interfacer.get_field(
    target_space::CC.Spaces.AbstractSpace,
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    quantity,
)
    return Interfacer.remap(
        target_space,
        Interfacer.get_field(sim, quantity),
        sim.remapping,
    )
end

"""
    get_oc_sim(sim)

Return the underlying `Oceananigans.Simulation` object for component models
that use Oceananigans under the hood.
"""
get_oc_sim(sim::OceananigansSimulation) = sim.ocean
get_oc_sim(sim::ClimaSeaIceSimulation) = sim.ice

"""
    Interfacer.sim_dt(sim::Union{OceananigansSimulation, ClimaSeaIceSimulation})

Return the simulation's timestep in seconds as a `Float64`.
"""
Interfacer.sim_dt(sim::Union{OceananigansSimulation, ClimaSeaIceSimulation}) =
    Float64(float(sim.model_Δt))

"""
    Interfacer.will_step(sim::Union{OceananigansSimulation, ClimaSeaIceSimulation}, t)

Return `true` if `Interfacer.step!(sim, t)` would take at least one step.
"""
function Interfacer.will_step(
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    t::Float64,
)
    oc_sim = get_oc_sim(sim)
    return (t - oc_sim.model.clock.time) >= Float64(sim.model_Δt)
end

function Interfacer.will_step(
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    t::ITime,
)
    oc_sim = get_oc_sim(sim)
    Δt_msec = date(t) - oc_sim.model.clock.time
    model_Δt_msec = counter(sim.model_Δt) * Dates.Millisecond(period(sim.model_Δt))
    return Δt_msec >= model_Δt_msec
end

"""
    Checkpointer.checkpoint_model_state(sim, comms_ctx, t, prev_checkpoint_t; output_dir)

Save the state of an Oceananigans-backed simulation to a JLD2 file at time `t`
(in seconds) using `Oceananigans.checkpoint`.

If a previous checkpoint exists, it is removed to avoid accumulating files.
A value of -1 for `prev_checkpoint_t` indicates there is no previous checkpoint.
"""
function Checkpointer.checkpoint_model_state(
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    comms_ctx::ClimaComms.AbstractCommsContext,
    t::Int,
    prev_checkpoint_t::Int;
    output_dir = "output",
)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving checkpoint $(nameof(sim)) model state to JLD2 on day $day second $sec"
    output_file = joinpath(output_dir, "checkpoint_$(nameof(sim))_$t.jld2")
    prev_checkpoint_file =
        joinpath(output_dir, "checkpoint_$(nameof(sim))_$(prev_checkpoint_t).jld2")
    Checkpointer.remove_checkpoint(prev_checkpoint_file, prev_checkpoint_t, comms_ctx)
    OC.checkpoint(get_oc_sim(sim); filepath = output_file)
    return nothing
end

"""
    Checkpointer.restart_model_state!(sim, input_file, comms_ctx)

Restore the state of an Oceananigans-backed simulation from a JLD2 checkpoint
file using `Oceananigans.set!`.

The coupler constructs `input_file` with a `.hdf5` extension; this method
replaces it with `.jld2` to match the format written by `checkpoint_model_state`.
"""
function Checkpointer.restart_model_state!(
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    input_file,
    comms_ctx,
)
    jld2_file = replace(input_file, ".hdf5" => ".jld2")
    ispath(jld2_file) || error("Oceananigans checkpoint file not found: $jld2_file")
    OC.set!(get_oc_sim(sim); checkpoint = jld2_file)
    return nothing
end

"""
    Checkpointer.restart_model_cache!(sim, input_file)

No-op for Oceananigans-backed simulations. All necessary state is restored via
`restart_model_state!`; there is no separate cache to restore.
"""
function Checkpointer.restart_model_cache!(
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    input_file,
)
    @warn "$(nameof(sim)) does not support restoring the model cache from a checkpoint. " *
          "The simulation cache will not be restored."
    return nothing
end
