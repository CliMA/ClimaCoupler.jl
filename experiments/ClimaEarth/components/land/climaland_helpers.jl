import ClimaCore as CC
import ClimaLand as CL
"""
    temp_anomaly_aquaplanet(coord)

Introduce a temperature IC anomaly for the aquaplanet case.
The values for this case follow the moist Held-Suarez setup of Thatcher &
Jablonowski (2016, eq. 6), consistent with ClimaAtmos aquaplanet.
"""
temp_anomaly_aquaplanet(coord) = 29 * exp(-coord.lat^2 / (2 * 26^2))

"""
    temp_anomaly_amip(coord)

Introduce a temperature IC anomaly for the AMIP case.
The values used in this case have been tuned to align with observed temperature
and result in stable simulations.
"""
temp_anomaly_amip(coord) = 40 * cosd(coord.lat)^4

"""
    make_land_domain(
        atmos_boundary_space::CC.Spaces.SpectralElementSpace2D,
        zlim::Tuple{FT, FT},
        nelements_vert::Int,) where {FT}

Creates the land model domain from the horizontal space of the atmosphere, and information
about the number of elements and extent of the vertical domain.
"""
function make_land_domain(
    atmos_boundary_space::CC.Spaces.SpectralElementSpace2D,
    zlim::Tuple{FT, FT},
    nelements_vert::Int,
) where {FT}
    @assert zlim[1] < zlim[2]
    depth = zlim[2] - zlim[1]

    mesh = CC.Spaces.topology(atmos_boundary_space).mesh

    radius = mesh.domain.radius
    nelements_horz = mesh.ne
    npolynomial = CC.Spaces.Quadratures.polynomial_degree(CC.Spaces.quadrature_style(atmos_boundary_space))
    nelements = (nelements_horz, nelements_vert)
    vertdomain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(FT(zlim[1])),
        CC.Geometry.ZPoint(FT(zlim[2]));
        boundary_names = (:bottom, :top),
    )

    vertmesh = CC.Meshes.IntervalMesh(vertdomain, CC.Meshes.Uniform(), nelems = nelements[2])
    verttopology = CC.Topologies.IntervalTopology(vertmesh)
    vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(verttopology)
    subsurface_space = CC.Spaces.ExtrudedFiniteDifferenceSpace(atmos_boundary_space, vert_center_space)
    subsurface_face_space = CC.Spaces.face_space(subsurface_space)
    space = (; surface = atmos_boundary_space, subsurface = subsurface_space, subsurface_face = subsurface_face_space)

    fields = CL.Domains.get_additional_coordinate_field_data(subsurface_space)

    if pkgversion(CL) < v"0.16.0"
        return CL.Domains.SphericalShell{FT}(radius, depth, nothing, nelements, npolynomial, space, fields)
    else
        return CL.Domains.SphericalShell{FT, typeof(space), typeof(fields)}(
            radius,
            depth,
            nothing,
            nelements,
            npolynomial,
            space,
            fields,
        )
    end
end
