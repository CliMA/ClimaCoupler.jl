"""
    to_node(pt::CC.Geometry.LatLongPoint)

Transform `LatLongPoint` into a tuple (long, lat, 0), where the 0 is needed because we only
care about the surface.
"""
@inline to_node(pt::CC.Geometry.LatLongPoint) = pt.long, pt.lat, zero(pt.lat)
# This next one is needed if we have "LevelGrid"
@inline to_node(pt::CC.Geometry.LatLongZPoint) = pt.long, pt.lat, zero(pt.lat)

"""
    map_interpolate(points, oc_field::OC.Field)

Interpolate the given 3D field onto the target points.

If the underlying grid does not contain a given point, return 0 instead.

Note: `map_interpolate` does not support interpolation from `Field`s defined on
`OrthogononalSphericalShellGrids` such as the `TripolarGrid`.

TODO: Use a non-allocating version of this function (simply replace `map` with `map!`)
"""
function map_interpolate(points, oc_field::OC.Field)
    loc = map(L -> L(), OC.Fields.location(oc_field))
    grid = oc_field.grid
    data = oc_field.data

    # TODO: There has to be a better way
    min_lat, max_lat = extrema(OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center()))

    map(points) do pt
        FT = eltype(pt)

        # The oceananigans grid does not cover the entire globe, so we should not
        # interpolate outside of its latitude bounds. Instead we return 0
        min_lat < pt.lat < max_lat || return FT(0)

        fᵢ = OC.Fields.interpolate(to_node(pt), data, loc, grid)
        convert(FT, fᵢ)::FT
    end
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

function Interfacer.remap(field::OC.Field, target_space)
    return map_interpolate(CC.Fields.coordinate_field(target_space), field)
end

function Interfacer.remap(operation::OC.AbstractOperations.AbstractOperation, target_space)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return Interfacer.remap(evaluated_field, target_space)
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
