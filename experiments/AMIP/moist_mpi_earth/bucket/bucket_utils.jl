"""
    get_land_temp(slab_sim::BucketSimulation)

Returns the surface temperature of the earth; 
a method for the bucket model 
when used as the land model.
"""
function get_land_temp(slab_sim::BucketSimulation)
    return slab_sim.integrator.p.bucket.T_sfc
end


"""
    get_land_temp_from_state(land_sim, u)

Returns the surface temperature of the earth, computed
from the state u.
"""
function get_land_temp_from_state(land_sim, u)
    face_space = ClimaLSM.Domains.obtain_face_space(land_sim.domain.domain.subsurface.space)
    N = ClimaCore.Spaces.nlevels(face_space)
    interp_c2f = ClimaCore.Operators.InterpolateC2F(
        top = ClimaCore.Operators.Extrapolate(),
        bottom = ClimaCore.Operators.Extrapolate(),
    )
    surface_space = land_sim.domain.domain.surface.space
    return ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(
            ClimaCore.Fields.level(interp_c2f.(u.bucket.T), ClimaCore.Utilities.PlusHalf(N - 1)),
        ),
        surface_space,
    )
end

"""
    get_land_roughness(slab_sim::BucketSimulation)

Returns the roughness length parameters of the bucket; 
a method for the bucket model 
when used as the land model.
"""
function get_land_roughness(slab_sim::BucketSimulation)
    return slab_sim.params.z_0m, slab_sim.params.z_0b
end

"""
   land_albedo(slab_sim::BucketSimulation)

Returns the surface albedo of the earth; 
a method for the bucket model 
when used as the land model.
"""
function land_albedo(slab_sim::BucketSimulation)
    α_land = surface_albedo.(slab_sim.integrator.p.bucket.α_sfc, slab_sim.params.albedo.α_snow, slab_sim.integrator.u.bucket.σS, slab_sim.params.σS_c)
    return parent(α_land)
end


"""
    get_land_q(slab_sim::Bucketimulation, _...)

Returns the surface specific humidity of the earth; 
a method for the bucket 
when used as the land model.
"""
function get_land_q(slab_sim::BucketSimulation, _...)
    return slab_sim.integrator.p.bucket.q_sfc
end

"""
    get_bucket_energy(bucket_sim)

Returns the volumetric internal energy of the bucket land model.
"""
function get_land_energy(bucket_sim::BucketSimulation, e_per_area)

    e_per_area .= zeros(axes(bucket_sim.integrator.u.bucket.W))
    soil_depth = FT = eltype(bucket_sim.integrator.u.bucket.W)
    ClimaCore.Fields.bycolumn(axes(bucket_sim.integrator.u.bucket.T)) do colidx
        e_per_area[colidx] .=
            bucket_sim.params.ρc_soil .* mean(bucket_sim.integrator.u.bucket.T[colidx]) .* bucket_sim.domain.soil_depth
    end

    e_per_area .+=
        -LSMP.LH_f0(bucket_sim.params.earth_param_set) .* LSMP.ρ_cloud_liq(bucket_sim.params.earth_param_set) .*
        bucket_sim.integrator.u.bucket.σS
    return e_per_area
end


"""
    make_lsm_domain(
        atmos_boundary_space::ClimaCore.Spaces.SpectralElementSpace2D,
        zlim::Tuple{FT, FT},
        nelements_vert::Int,) where {FT}

Creates the LSM Domain from the horizontal space of the atmosphere, and information
about the number of elements and extent of the vertical domain.
"""
function make_lsm_domain(
    atmos_boundary_space::ClimaCore.Spaces.SpectralElementSpace2D,
    zlim::Tuple{FT, FT},
    nelements_vert::Int,
) where {FT}
    @assert zlim[1] < zlim[2]
    height = zlim[2] - zlim[1]
    radius = atmos_boundary_space.topology.mesh.domain.radius
    npolynomial = ClimaCore.Spaces.Quadratures.polynomial_degree(atmos_boundary_space.quadrature_style)
    nelements_horz = atmos_boundary_space.topology.mesh.ne
    nelements = (nelements_horz, nelements_vert)
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(zlim[1])),
        ClimaCore.Geometry.ZPoint(FT(zlim[2]));
        boundary_tags = (:bottom, :top),
    )

    vertmesh = ClimaCore.Meshes.IntervalMesh(vertdomain, ClimaCore.Meshes.Uniform(), nelems = nelements[2])
    verttopology = ClimaCore.Topologies.IntervalTopology(vertmesh)
    vert_center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(verttopology)
    subsurface_space = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(atmos_boundary_space, vert_center_space)

    surface_domain = ClimaLSM.Domains.SphericalSurface{FT, typeof(atmos_boundary_space)}(
        radius,
        nelements[1],
        npolynomial,
        atmos_boundary_space,
    )
    subsurface_domain = ClimaLSM.Domains.SphericalShell{FT, typeof(subsurface_space)}(
        radius,
        height,
        nelements,
        npolynomial,
        subsurface_space,
    )
    return ClimaLSM.Domains.LSMSphericalShellDomain{FT, typeof(subsurface_space), typeof(atmos_boundary_space)}(
        subsurface_domain,
        surface_domain,
    )
end
