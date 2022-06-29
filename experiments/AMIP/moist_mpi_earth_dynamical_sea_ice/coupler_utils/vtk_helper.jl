## Plot VTK output on the 2D cubed-sphere, building on from ClimaCoreWorkshops
using ClimaCoreVTK

# Plot the time series

"""
    convert_field_to_type(space_old, space_new, field)
- used specifically to convert a Field from Float32 to Float64, which is acceptable for ClimaCoreVTK.writepvd
"""
function convert_field_to_type(space_new, field_old)
    fieldT = ones(space_new)
    fieldT_data = convert.(eltype(space_new), parent(field_old))
    parent(fieldT) .= fieldT_data
    return fieldT
end

"""
    retype_space(space::ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace{CellCenter}, T)
- used specifically to convert a space from Float32 to Float64, which feed into `convert_field_to_type`
"""
function retype_space(space::ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace{ClimaCore.Spaces.CellCenter}, T)
    R = space.horizontal_space.topology.mesh.domain.radius
    ne = space.horizontal_space.topology.mesh.ne
    Nq = Spaces.Quadratures.polynomial_degree(space.horizontal_space.quadrature_style) + 1

    horizontal_meshT = cubed_sphere_mesh(; radius = T(R), h_elem = ne)
    h_spaceT = make_horizontal_space(horizontal_meshT, quad, nothing)

    z_stretch = Meshes.GeneralizedExponentialStretching(T(500), T(5000)) # TODO
    z_max = space.vertical_topology.mesh.domain.coord_max.z
    z_elem = 10  # TODO
    spaceT, _ = make_hybrid_spaces(h_spaceT, T(z_max), z_elem, z_stretch)
    return spaceT
end

function retype_space(space::ClimaCore.Spaces.SpectralElementSpace2D, T)
    R = space.topology.mesh.domain.radius
    ne = space.topology.mesh.ne
    Nq = Spaces.Quadratures.polynomial_degree(space.quadrature_style) + 1

    horizontal_meshT = cubed_sphere_mesh(; radius = T(R), h_elem = ne)
    h_spaceT = make_horizontal_space(horizontal_meshT, quad, nothing)

    return h_spaceT
end
