"""
    get_land_temp_from_state(land_sim, u)

Returns the surface temperature of the earth, computed
from the state u.
"""
function get_land_temp_from_state(land_sim, u)
    return ClimaLSM.surface_temperature(land_sim.model, u, land_sim.integrator.p, land_sim.integrator.t)
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
            bucket_sim.model.parameters.ρc_soil .* mean(bucket_sim.integrator.u.bucket.T[colidx]) .*
            bucket_sim.domain.soil_depth
    end

    e_per_area .+=
        -LSMP.LH_f0(bucket_sim.model.parameters.earth_param_set) .*
        LSMP.ρ_cloud_liq(bucket_sim.model.parameters.earth_param_set) .* bucket_sim.integrator.u.bucket.σS
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

function update_calculated_fluxes_point!(sim::BucketSimulation, fields, colidx)
    (; ρ_sfc, F_shf, F_lhf, F_evap) = fields
    ρ_liq = (LSMP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    @. sim.integrator.p.bucket.ρ_sfc[colidx] = ρ_sfc
    @. sim.integrator.p.bucket.turbulent_energy_flux[colidx] = F_shf + F_lhf
    @. sim.integrator.p.bucket.evaporation[colidx] = F_evap / ρ_liq
    return nothing
end

function update!(sim::BucketSimulation, ::Val{:net_radiation}, field)
    @. sim.integrator.p.bucket.R_n .= field
end
function update!(sim::BucketSimulation, ::Val{:precipitation_liquid}, field)
    @. sim.integrator.p.bucket.P_liq .= field
end
function update!(sim::BucketSimulation, ::Val{:precipitation_snow}, field)
    @. sim.integrator.p.bucket.P_snow .= field
end

function update!(sim::BucketSimulation, ::Val{:ρ_sfc}, field)
    parent(land_sim.integrator.p.bucket.ρ_sfc) .= parent(field)
end

#  combining means atmos needs to know about the coupler interface.
# function update_sim!(sim::BucketSimulation, csf, fraction)

#     @. sim.integrator.p.bucket.R_n .= csf.F_R_sfc * fraction
#     @. sim.integrator.p.bucket.P_liq .= csf.P_liq * fraction
#     @. sim.integrator.p.bucket.P_snow .= csf.P_snow * fraction
#     @. land_sim.integrator.p.bucket.ρ_sfc .= csf.ρ_sfc * fraction
# end

# parent(land_sim.integrator.p.bucket.ρ_sfc) .= parent(csf.ρ_sfc)
# parent(land_sim.integrator.p.bucket.turbulent_energy_flux) .=
#     parent(Regridder.binary_mask.(land_fraction, threshold = eps()) .* csf.F_A)
# ρ_liq = (LSMP.ρ_cloud_liq(land_sim.model.parameters.earth_param_set))
# parent(land_sim.integrator.p.bucket.evaporation) .=
#     parent(Regridder.binary_mask.(land_fraction, threshold = eps()) .* csf.F_evap ./ ρ_liq)
# parent(land_sim.integrator.p.bucket.R_n) .=
#     parent(Regridder.binary_mask.(land_fraction, threshold = eps()) .* csf.F_R_sfc)
# parent(land_sim.integrator.p.bucket.P_liq) .= FT(-1.0) .* parent(csf.P_liq) # land expects this to be positive
# parent(land_sim.integrator.p.bucket.P_snow) .= FT(0.0) .* parent(csf.P_snow)


"""
    get_...(sim::BucketSimulation, colidx)

Returns the surface the respective variables of the earth at column index `colidx`;
a method for the bucket when used as the land model.
"""

get_temperature(sim::BucketSimulation) = ClimaLSM.surface_temperature( sim.model, sim.integrator.u, sim.integrator.p, sim.integrator.t)
get_humidity(sim::BucketSimulation) = ClimaLSM.surface_specific_humidity(sim.model, sim.integrator.u, sim.integrator.p)
get_z0m(sim::BucketSimulation) = sim.model.parameters.z_0m
get_z0b(sim::BucketSimulation) = sim.model.parameters.z_0b
get_beta(sim::BucketSimulation) = ClimaLSM.surface_evaporative_scaling(sim.model, sim.integrator.u, sim.integrator.p)
get_albedo(sim::BucketSimulation) = ClimaLSM.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
get_area_fraction(sim::BucketSimulation) = sim.area_fraction

get_temperature_point(sim::BucketSimulation, colidx) = get_temperature(sim)[colidx]
get_humidity_point(sim::BucketSimulation, colidx) = get_humidity(sim)[colidx]
get_z0m_point(sim::BucketSimulation, colidx) = get_z0m(sim)
get_z0b_point(sim::BucketSimulation, colidx) = get_z0b(sim)
get_beta_point(sim::BucketSimulation, colidx) = get_beta(sim)[colidx]
get_albedo_point(sim::BucketSimulation, colidx) = get_albedo(sim)[colidx]
get_heat_transfer_coefficient_point(sim::BucketSimulation, colidx) = sim.integrator.p.Ch
get_drag_transfer_coefficient_point(sim::BucketSimulation, colidx) = sim.integrator.p.Cd

issaturated(::BucketSimulation, q) = isnan(parent(q)[1])

reinit!(sim::BucketSimulation) = reinit!(sim.integrator)