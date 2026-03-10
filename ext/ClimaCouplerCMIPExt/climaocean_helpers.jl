"""
    to_node(pt::CC.Geometry.LatLongPoint)

Transform `LatLongPoint` into a tuple (long, lat, 0), where the 0 is needed because we only
care about the surface.
"""
@inline to_node(pt::CC.Geometry.LatLongPoint) = pt.long, pt.lat, zero(pt.lat)
# This next one is needed if we have "LevelGrid"
@inline to_node(pt::CC.Geometry.LatLongZPoint) = pt.long, pt.lat, zero(pt.lat)

"""
    map_interpolate!(target_field::CC.Fields.Field, source_field::OC.Field)

Interpolate the given 3D field onto the target field, modifying the target field in place.

If the underlying grid does not contain a given point, writes 0 instead.

Note: `map_interpolate!` does not support interpolation from `Field`s defined on
`OrthogononalSphericalShellGrids` such as the `TripolarGrid`.
"""
function map_interpolate!(target_field::CC.Fields.Field, source_field::OC.Field)
    points = CC.Fields.coordinate_field(axes(target_field))
    loc = map(L -> L(), OC.Fields.location(source_field))
    grid = source_field.grid
    data = source_field.data

    # TODO: There has to be a better way
    min_lat, max_lat = extrema(OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center()))

    # Use map! on the field directly, which handles complex data layouts correctly
    map!(target_field, points) do pt
        FT = eltype(pt)

        # The oceananigans grid does not cover the entire globe, so we should not
        # interpolate outside of its latitude bounds. Instead we return 0
        min_lat < pt.lat < max_lat || return FT(0)

        fᵢ = OC.Fields.interpolate(to_node(pt), data, loc, grid)
        return convert(FT, fᵢ)::FT
    end
    return nothing
end

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

function Interfacer.remap!(target_field::CC.Fields.Field, source_field::OC.Field)
    map_interpolate!(target_field, source_field)
    return nothing
end

function Interfacer.remap!(
    target_field::CC.Fields.Field,
    operation::OC.AbstractOperations.AbstractOperation,
)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    Interfacer.remap!(target_field, evaluated_field)
    return nothing
end

function Interfacer.remap(target_space, field::OC.Field)
    # Allocate target field and call remap!
    target_field = CC.Fields.zeros(target_space)
    Interfacer.remap!(target_field, field)
    return target_field
end

function Interfacer.remap(target_space, operation::OC.AbstractOperations.AbstractOperation)
    # Allocate target field and call remap!
    target_field = CC.Fields.zeros(target_space)
    Interfacer.remap!(target_field, operation)
    return target_field
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
