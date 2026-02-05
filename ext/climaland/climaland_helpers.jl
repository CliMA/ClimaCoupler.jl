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
        depth::FT,
        toml_dict::CP.ParamDict;
        nelements::Tuple{Int, Int} = (101, 15),
        dz_tuple::Tuple{FT, FT} = FT.((10.0, 0.05)),
    ) where {FT}
Creates a land model domain with the given number of vertical elements
and vertical depth.
"""
function make_land_domain(
    depth::FT,
    toml_dict::CP.ParamDict;
    nelements::Tuple{Int, Int} = (101, 15),
    dz_tuple::Tuple{FT, FT} = FT.((10.0, 0.05)),
) where {FT}
    radius = toml_dict["planet_radius"] # in meters
    npolynomial = 0
    domain = CL.Domains.SphericalShell(; radius, depth, nelements, npolynomial, dz_tuple)
    return domain
end

"""
     make_land_domain(
         shared_surface_space::CC.Spaces.SpectralElementSpace2D,
         depth::FT;
         nelements_vert::Int = 15,
         dz_tuple::Tuple{FT, FT} = FT.((10.0, 0.05)),
         ) where {FT}

 Creates the land model domain from the input shared surface space and information
 about the number of elements and extent of the vertical domain.
 """
function make_land_domain(
    shared_surface_space::CC.Spaces.SpectralElementSpace2D,
    depth::FT;
    nelements_vert::Int = 15,
    dz_tuple::Tuple{FT, FT} = FT.((10.0, 0.05)),
) where {FT}
    mesh = CC.Spaces.topology(shared_surface_space).mesh

    radius = mesh.domain.radius
    nelements_horz = mesh.ne
    npolynomial = CC.Spaces.Quadratures.polynomial_degree(
        CC.Spaces.quadrature_style(shared_surface_space),
    )
    nelements = (nelements_horz, nelements_vert)
    vertdomain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(FT(-depth)),
        CC.Geometry.ZPoint(FT(0));
        boundary_names = (:bottom, :top),
    )

    vertmesh = CC.Meshes.IntervalMesh(
        vertdomain,
        CC.Meshes.GeneralizedExponentialStretching{FT}(dz_tuple[1], dz_tuple[2]);
        nelems = nelements_vert,
        reverse_mode = true,
    )
    verttopology = CC.Topologies.IntervalTopology(vertmesh)
    vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(verttopology)
    subsurface_space =
        CC.Spaces.ExtrudedFiniteDifferenceSpace(shared_surface_space, vert_center_space)
    subsurface_face_space = CC.Spaces.face_space(subsurface_space)
    space = (;
        surface = shared_surface_space,
        subsurface = subsurface_space,
        subsurface_face = subsurface_face_space,
    )

    fields = CL.Domains.get_additional_coordinate_field_data(subsurface_space)

    return CL.Domains.SphericalShell{FT, typeof(space), typeof(fields)}(
        radius,
        depth,
        dz_tuple,
        nelements,
        npolynomial,
        space,
        fields,
    )
end
