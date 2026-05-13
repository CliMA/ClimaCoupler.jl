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


### Helper functions for non-conservative remapping
"""
    to_node(pt::CC.Geometry.LatLongPoint)

Transform `LatLongPoint` into a tuple `(long, lat, 0)`, where the trailing 0 is
needed because we only care about the surface. Helper used by `map_interpolate!`
for non-conservative OC → CC remapping on `LatitudeLongitudeGrid`s.
"""
@inline to_node(pt::CC.Geometry.LatLongPoint) = pt.long, pt.lat, zero(pt.lat)
# This next one is needed if we have "LevelGrid"
@inline to_node(pt::CC.Geometry.LatLongZPoint) = pt.long, pt.lat, zero(pt.lat)

"""
    map_interpolate!(target_field::CC.Fields.Field, source_field::OC.Field)

Interpolate the given 3D Oceananigans field onto the ClimaCore target field,
modifying the target field in place. If the underlying grid does not contain a
given point, writes 0 instead.

This is the non-conservative OC → CC path used with `LatitudeLongitudeGrid`s; it
does **not** support `OrthogonalSphericalShellGrid`s like `TripolarGrid`.
"""
function map_interpolate!(target_field::CC.Fields.Field, source_field::OC.Field)
    points = CC.Fields.coordinate_field(axes(target_field))
    loc = map(L -> L(), OC.Fields.location(source_field))
    grid = source_field.grid
    data = source_field.data

    # TODO: There has to be a better way to compute the latitude bounds.
    min_lat, max_lat = extrema(OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center()))

    # Use map! on the field directly, which handles complex data layouts correctly
    map!(target_field, points) do pt
        FT = eltype(pt)

        # The oceananigans grid does not cover the entire globe, so we should
        # not interpolate outside of its latitude bounds; instead, return 0.
        min_lat < pt.lat < max_lat || return FT(0)

        fᵢ = OC.Fields.interpolate(to_node(pt), data, loc, grid)
        return convert(FT, fᵢ)::FT
    end
    return nothing
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

Convert the covariant vector components `ρτxz` and `ρτyz` from the
contravariant basis (as they are output by the surface flux calculation)
to the Cartesian basis. These are now in an extrinsic coordinate system
that can be rotated onto the ocean/sea ice grid by `_rotate_vector!`.
"""
function contravariant_to_cartesian!(ρτ_flux_uv, ρτxz, ρτyz)
    # Get the local geometry of the boundary space
    local_geometry = CC.Fields.local_geometry_field(ρτxz)

    # Get the vector components in the CT1 and CT2 directions
    xz = @. CT12(CT1(unit_basis_vector_data(CT1, local_geometry)), local_geometry)
    yz = @. CT12(CT2(unit_basis_vector_data(CT2, local_geometry)), local_geometry)

    # Convert the vector components to a UVVector on the Cartesian basis
    @. ρτ_flux_uv = CC.Geometry.UVVector(ρτxz * xz + ρτyz * yz, local_geometry)
    return nothing
end

# Define shorthands for ClimaCore types
const CT1 = CC.Geometry.Contravariant1Vector
const CT2 = CC.Geometry.Contravariant2Vector
const CT12 = CC.Geometry.Contravariant12Vector

"""
    unit_basis_vector_data(type, local_geometry)

The component of the vector of the specified type with length 1 in physical units.
The type should correspond to a vector with only one component, i.e., a basis vector.

Helper function used only in `contravariant_to_cartesian!`.
"""
function unit_basis_vector_data(::Type{V}, local_geometry) where {V}
    FT = CC.Geometry.undertype(typeof(local_geometry))
    return FT(1) / CC.Geometry._norm(V(FT(1)), local_geometry)
end

### Extensions of Interfacer.jl remapping functions for Oceananigans fields/grids
# Non-allocating ClimaCore -> Oceananigans remap
function Interfacer.remap!(
    target_field::OC.Field,
    source_field::CC.Fields.Field,
    remapping,
)
    return _remap_helper!(
        target_field,
        source_field,
        remapping,
        remapping.remapper_cc_to_oc,
    )
end

# Non-allocating Oceananigans -> ClimaCore remap
function Interfacer.remap!(
    target_field::CC.Fields.Field,
    source_field::OC.Field,
    remapping,
)
    return _remap_helper!(
        target_field,
        source_field,
        remapping,
        remapping.remapper_cc_to_oc,
    )
end

# CC -> OC, conservative
function _remap_helper!(
    target_field::OC.Field,
    source_field::CC.Fields.Field,
    remapping,
    ::CR.Regridder,
)
    get_ConservativeRegriddingCCExt().get_value_per_element!(
        remapping.value_per_element_cc,
        source_field,
        remapping.field_ones_cc,
    )

    # Regrid into a contiguous (Nx, Ny) staging buffer, then copy into the OC.Field's
    # (strided, non-contiguous) interior.
    CR.regrid!(
        vec(remapping.remap_scratch_arr),
        remapping.remapper_cc_to_oc,
        remapping.value_per_element_cc,
    )
    # Get the index of the top level (surface); 1 for 2D fields, Nz for 3D fields
    z = size(target_field, 3)
    OC.interior(target_field, :, :, z) .= remapping.remap_scratch_arr
    return nothing
end

# CC -> OC, non-conservative
function _remap_helper!(
    target_field::OC.Field,
    source_field::CC.Fields.Field,
    remapping,
    remapper_cc_to_oc::CC.Remapping.Remapper,
)
    # `CC.Remapping.interpolate!` is written to fill a plain (contiguous) 2D buffer,
    # so we route through `remap_scratch_arr` and then copy into the OC.Field's
    # interior, rather than handing it a strided `OC.interior` view directly.
    CC.Remapping.interpolate!(
        remapping.remap_scratch_arr,
        remapper_cc_to_oc,
        source_field,
    )
    # Get the index of the top level (surface); 1 for 2D fields, Nz for 3D fields
    z = size(target_field, 3)
    OC.interior(target_field, :, :, z) .= remapping.remap_scratch_arr
    return nothing
end

# OC -> CC, conservative
function _remap_helper!(
    target_field::CC.Fields.Field,
    source_field::OC.Field,
    remapping,
    ::CR.Regridder,
)
    # Stage the OC.Field's strided interior into a contiguous (Nx, Ny) buffer
    # before regridding. This guarantees `CR.regrid!` receives a plain
    # Vector/CuVector view (via `vec`), rather than a `ReshapedArray` wrapper over
    # the halo'd parent's strided SubArray.
    z = size(source_field, 3)
    remapping.remap_scratch_arr .= OC.interior(source_field, :, :, z)

    # Regrid the source field; results land in the per-element CC vector. The
    # regridder is built CC → OC, so we use `transpose` to go the other way.
    CR.regrid!(
        remapping.value_per_element_cc,
        transpose(remapping.remapper_cc_to_oc),
        vec(remapping.remap_scratch_arr),
    )

    # Convert the vector of remapped values to a ClimaCore Field with one value per element
    get_ConservativeRegriddingCCExt().set_value_per_element!(
        target_field,
        remapping.value_per_element_cc,
    )
    return nothing
end

# OC -> CC, non-conservative
function _remap_helper!(
    target_field::CC.Fields.Field,
    source_field::OC.Field,
    remapping,
    ::CC.Remapping.Remapper,
)
    map_interpolate!(target_field, source_field)
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

# Handle the case of remapping a scalar number to a ClimaCore space
Interfacer.remap!(target_field::CC.Fields.Field, source_field::Number, remapping) =
    Interfacer.remap!(target_field, source_field)
Interfacer.remap(target_space::CC.Spaces.AbstractSpace, source_num::Number, remapping) =
    Interfacer.remap(target_space, source_num, remapping)

# Handle the case of remapping the area fraction field, which is a ClimaCore Field
Interfacer.remap!(target_field::CC.Fields.Field, source_field::CC.Fields.Field, remapping) =
    Interfacer.remap!(target_field, source_field)

# Extend Interfacer.get_field to allow automatic remapping to the target space
function Interfacer.get_field!(
    target_field,
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    quantity,
)
    Interfacer.remap!(target_field, Interfacer.get_field(sim, quantity), sim.remapping)
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
