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

"""
    contravariant_to_cartesian!(ρτ_flux_uv, ρτxz, ρτyz)

Convert contravariant momentum-flux components `ρτxz` and `ρτyz` (CT1/CT2
scalars on the boundary space, as returned by `SurfaceFluxes`) to extrinsic
`UVVector` components using the local geometry at each node.
"""
function contravariant_to_cartesian!(ρτ_flux_uv, ρτxz, ρτyz)
    local_geometry = CC.Fields.local_geometry_field(ρτxz)
    xz = @. CT12(CT1(unit_basis_vector_data(CT1, local_geometry)), local_geometry)
    yz = @. CT12(CT2(unit_basis_vector_data(CT2, local_geometry)), local_geometry)
    @. ρτ_flux_uv = CC.Geometry.UVVector(ρτxz * xz + ρτyz * yz, local_geometry)
    return nothing
end

const CT1 = CC.Geometry.Contravariant1Vector
const CT2 = CC.Geometry.Contravariant2Vector
const CT12 = CC.Geometry.Contravariant12Vector

function unit_basis_vector_data(::Type{V}, local_geometry) where {V}
    FT = CC.Geometry.undertype(typeof(local_geometry))
    return FT(1) / CC.Geometry._norm(V(FT(1)), local_geometry)
end

"""
    remap_contravariant_momentum_fluxes_to_oc!(scratch_oc1, scratch_oc2, temp_uv_vec,
                                               remapping, F_turb_ρτxz, F_turb_ρτyz)

Convert contravariant momentum fluxes on the boundary space to `UVVector`,
then conservatively remap the extrinsic components onto OC scratch fields.

!!! warning "Per-CC-element constant output"
    This function routes through `CR.regrid!` (SE→FV principled operator) on a
    field that was built by broadcasting a per-element scalar average across all
    Nq² GLL nodes.  The result is constant within each CC element on the OC grid.
    Prefer [`scatter_contravariant_momentum_to_oc!`](@ref) with precomputed polygon
    geometry factors to obtain per-OC-cell resolution.
"""
function remap_contravariant_momentum_fluxes_to_oc!(
    scratch_oc1,
    scratch_oc2,
    temp_uv_vec,
    remapping,
    F_turb_ρτxz,
    F_turb_ρτyz,
)
    contravariant_to_cartesian!(temp_uv_vec, F_turb_ρτxz, F_turb_ρτyz)
    Interfacer.remap!(scratch_oc1, temp_uv_vec.components.data.:1, remapping)
    Interfacer.remap!(scratch_oc2, temp_uv_vec.components.data.:2, remapping)
    return nothing
end

"""
    precompute_intersection_momentum_geometry(ig, boundary_space_cpu)

Precompute per-intersection-polygon UV (east-north) components of the unit
contravariant CT1 and CT2 stress vectors.

The result is a `NamedTuple` with four per-polygon `Vector{FT}` fields:
- `ct1_u`, `ct1_v`: U and V components arising from unit CT1 stress
- `ct2_u`, `ct2_v`: U and V components arising from unit CT2 stress

These are used by [`scatter_contravariant_momentum_to_oc!`](@ref) to convert
per-polygon contravariant stress directly to Cartesian UV at each intersection
polygon, avoiding the CC-element aggregation that would otherwise collapse the
OC-resolution information.

The geometry factors are static (they depend only on the mesh, not the state)
and should be precomputed once during `construct_remapper` setup.
"""
function precompute_intersection_momentum_geometry(ig::IntersectionGrid, boundary_space_cpu)
    FT = CC.Spaces.undertype(boundary_space_cpu)
    CRExt = get_ConservativeRegriddingCCExt()

    # Scratch fields for contravariant→UV evaluation.
    tmp1   = CC.Fields.zeros(boundary_space_cpu)
    tmp2   = CC.Fields.zeros(boundary_space_cpu)
    tmp_uv = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space_cpu)

    n_int = ig.n_intersections
    ct1_u = zeros(FT, n_int)
    ct1_v = zeros(FT, n_int)
    ct2_u = zeros(FT, n_int)
    ct2_v = zeros(FT, n_int)

    tmp1 .= one(FT)
    tmp2 .= zero(FT)
    contravariant_to_cartesian!(tmp_uv, tmp1, tmp2)
    gather_cc_nodal_to_intersection!(ct1_u, ig, copy(CRExt.se_field_to_vec(tmp_uv.components.data.:1)))
    gather_cc_nodal_to_intersection!(ct1_v, ig, copy(CRExt.se_field_to_vec(tmp_uv.components.data.:2)))

    # Unit CT2 stress (0,1) → UV at every GLL node.
    tmp1 .= zero(FT)
    tmp2 .= one(FT)
    contravariant_to_cartesian!(tmp_uv, tmp1, tmp2)
    gather_cc_nodal_to_intersection!(ct2_u, ig, copy(CRExt.se_field_to_vec(tmp_uv.components.data.:1)))
    gather_cc_nodal_to_intersection!(ct2_v, ig, copy(CRExt.se_field_to_vec(tmp_uv.components.data.:2)))

    return (; ct1_u, ct1_v, ct2_u, ct2_v)
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

# Non-allocating Oceananigans operation -> ClimaCore remap
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

# Note: the `Interfacer.get_field!` / `Interfacer.get_field` methods that
# remap onto a target space for `Union{OceananigansSimulation,
# ClimaSeaIceSimulation}` live in `clima_seaice.jl`, which is included after
# the type definitions of both simulations.

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
