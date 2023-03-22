"""
    get_land_temp(slab_sim::BucketSimulation)

Returns the surface temperature of the earth;
a method for the bucket model
when used as the land model.
"""
function get_land_temp(slab_sim::BucketSimulation)
    return ClimaLSM.surface_temperature(
        slab_sim.model,
        slab_sim.integrator.u,
        slab_sim.integrator.p,
        slab_sim.integrator.t,
    )
end


"""
    get_land_temp_from_state(land_sim, u)

Returns the surface temperature of the earth, computed
from the state u.
"""
function get_land_temp_from_state(land_sim, u)
    return ClimaLSM.surface_temperature(land_sim.model, u, land_sim.integrator.p, land_sim.integrator.t)
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
    return ClimaLSM.surface_albedo(slab_sim.model, slab_sim.integrator.u, slab_sim.integrator.p)
end

"""
   land_beta(slab_sim::BucketSimulation)

Returns the surface evaporative scaling factor over land;
a method for the bucket model when used as the land model.
Note that this is slightly different from the coupler's β,
which includes the scaling factor over non-land surfaces.
"""
function land_beta(slab_sim::BucketSimulation)
    return ClimaLSM.surface_evaporative_scaling(slab_sim.model, slab_sim.integrator.u, slab_sim.integrator.p)
end


"""
    get_land_q(slab_sim::Bucketimulation, _...)

Returns the surface specific humidity of the earth;
a method for the bucket
when used as the land model.
"""
function get_land_q(slab_sim::BucketSimulation, _...)
    return ClimaLSM.surface_specific_humidity(slab_sim.model, slab_sim.integrator.u, slab_sim.integrator.p)
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
