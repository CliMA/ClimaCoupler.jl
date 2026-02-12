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
    @. ρτ_flux_uv = @. CC.Geometry.UVVector(ρτxz * xz + ρτyz * yz, local_geometry)
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
function Interfacer.remap!(target_field::OC.Field, source_field::CC.Fields.Field, remapping)
    ConservativeRegriddingCCExt.get_value_per_element!(
        remapping.value_per_element_cc,
        source_field,
        remapping.field_ones_cc,
    )

    # Get the index of the top level (surface); 1 for 2D fields, Nz for 3D fields
    z = size(target_field, 3)
    dst = vec(OC.interior(target_field, :, :, z))
    src = remapping.value_per_element_cc

    # Regrid the source field to the target field
    CR.regrid!(dst, transpose(remapping.remapper_oc_to_cc), src)
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

# Non-allocating Oceananigans Field -> ClimaCore remap
function Interfacer.remap!(target_field::CC.Fields.Field, source_field::OC.Field, remapping)
    # Get the index of the top level (surface); 1 for 2D fields, Nz for 3D fields
    z = size(source_field, 3)
    src = vec(OC.interior(source_field, :, :, z))

    # Store the remapped FV values in a vector of length equal to the number of elements in the target space
    dst = remapping.value_per_element_cc

    # Regrid the source field to the target field
    CR.regrid!(dst, remapping.remapper_oc_to_cc, src)

    # Convert the vector of remapped values to a ClimaCore Field with one value per element
    ConservativeRegriddingCCExt.set_value_per_element!(target_field, dst)
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
    sim::OceananigansSimulation,
    quantity,
)
    return Interfacer.remap(
        target_space,
        Interfacer.get_field(sim, quantity),
        sim.remapping,
    )
end
