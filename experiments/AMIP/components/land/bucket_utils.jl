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
    depth = zlim[2] - zlim[1]

    radius = ClimaCore.Spaces.topology(atmos_boundary_space).mesh.domain.radius
    nelements_horz = ClimaCore.Spaces.topology(atmos_boundary_space).mesh.ne
    npolynomial =
        ClimaCore.Spaces.Quadratures.polynomial_degree(ClimaCore.Spaces.quadrature_style(atmos_boundary_space))
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
    space = (; surface = atmos_boundary_space, subsurface = subsurface_space)

    return ClimaLSM.Domains.SphericalShell{FT}(radius, depth, nothing, nelements, npolynomial, space)
end

# extensions required by Interfacer
get_field(sim::BucketSimulation, ::Val{:surface_temperature}) =
    ClimaLSM.surface_temperature(sim.model, sim.integrator.u, sim.integrator.p, sim.integrator.t)
get_field(sim::BucketSimulation, ::Val{:surface_humidity}) =
    ClimaLSM.surface_specific_humidity(sim.model, sim.integrator.u, sim.integrator.p, sim.integrator.t)
get_field(sim::BucketSimulation, ::Val{:roughness_momentum}) = sim.model.parameters.z_0m
get_field(sim::BucketSimulation, ::Val{:roughness_buoyancy}) = sim.model.parameters.z_0b
get_field(sim::BucketSimulation, ::Val{:beta}) =
    ClimaLSM.surface_evaporative_scaling(sim.model, sim.integrator.u, sim.integrator.p)
get_field(sim::BucketSimulation, ::Val{:albedo}) =
    ClimaLSM.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
get_field(sim::BucketSimulation, ::Val{:area_fraction}) = sim.area_fraction
get_field(sim::BucketSimulation, ::Val{:air_density}) = sim.integrator.p.bucket.ρ_sfc

function update_field!(sim::BucketSimulation, ::Val{:turbulent_energy_flux}, field)
    parent(sim.integrator.p.bucket.turbulent_energy_flux) .= parent(field)
end
function update_field!(sim::BucketSimulation, ::Val{:turbulent_moisture_flux}, field)
    ρ_liq = (LSMP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    parent(sim.integrator.p.bucket.evaporation) .= parent(field ./ ρ_liq)
end
function update_field!(sim::BucketSimulation, ::Val{:radiative_energy_flux}, field)
    parent(sim.integrator.p.bucket.R_n) .= parent(field)
end
function update_field!(sim::BucketSimulation, ::Val{:liquid_precipitation}, field)
    ρ_liq = (LSMP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    parent(sim.integrator.p.drivers.P_liq) .= parent(field ./ ρ_liq)
end
function update_field!(sim::BucketSimulation, ::Val{:snow_precipitation}, field)
    ρ_liq = (LSMP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    parent(sim.integrator.p.drivers.P_snow) .= parent(field ./ ρ_liq)
end

function update_field!(sim::BucketSimulation, ::Val{:air_density}, field)
    parent(sim.integrator.p.bucket.ρ_sfc) .= parent(field)
end

step!(sim::BucketSimulation, t) = step!(sim.integrator, t - sim.integrator.t, true)

reinit!(sim::BucketSimulation) = reinit!(sim.integrator)

# extensions required by FluxCalculator (partitioned fluxes)
function update_turbulent_fluxes_point!(sim::BucketSimulation, fields::NamedTuple, colidx::Fields.ColumnIndex)
    (; F_turb_energy, F_turb_moisture) = fields
    sim.integrator.p.bucket.turbulent_energy_flux[colidx] .= F_turb_energy
    sim.integrator.p.bucket.evaporation[colidx] .=
        F_turb_moisture ./ LSMP.ρ_cloud_liq(sim.model.parameters.earth_param_set)
    return nothing
end

# extension of FluxCalculator.surface_thermo_state, overriding the saturated-surface default
function surface_thermo_state(
    sim::BucketSimulation,
    thermo_params::TD.Parameters.ThermodynamicsParameters,
    thermo_state_int,
    colidx::ClimaCore.Fields.ColumnIndex,
)

    T_sfc = get_field(sim, Val(:surface_temperature), colidx)
    # Note that the surface air density, ρ_sfc, is computed using the atmospheric state at the first level and making ideal gas
    # and hydrostatic balance assumptions. The land model does not compute the surface air density so this is
    # a reasonable stand-in.
    ρ_sfc = get_field(sim, Val(:air_density), colidx)
    q_sfc = get_field(sim, Val(:surface_humidity), colidx) # already calculated in rhs! (cache)
    @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
end

"""
    get_model_state_vector(sim::BucketSimulation)

Extension of Checkpointer.get_model_state_vector to get the model state.
"""
function get_model_state_vector(sim::BucketSimulation)
    return sim.integrator.u.bucket
end

"""
    get_field(bucket_sim::BucketSimulation, ::Val{:energy})

Extension of Interfacer.get_field that provides the total energy contained in the bucket, including the latent heat due to snow melt.
"""
function get_field(bucket_sim::BucketSimulation, ::Val{:energy})
    # required by ConservationChecker
    e_per_area = zeros(axes(bucket_sim.integrator.u.bucket.W))
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
    get_field(bucket_sim::BucketSimulation, ::Val{:water})

Extension of Interfacer.get_field that provides the total water contained in the bucket, including the liquid water in snow.
"""
function get_field(bucket_sim::BucketSimulation, ::Val{:water})
    ρ_cloud_liq = ClimaLSM.LSMP.ρ_cloud_liq(bucket_sim.model.parameters.earth_param_set)
    return
    @. (bucket_sim.integrator.u.bucket.σS + bucket_sim.integrator.u.bucket.W + bucket_sim.integrator.u.bucket.Ws) *
       ρ_cloud_liq  # kg water / m2
end

"""
    get_land_temp_from_state(land_sim, u)
Returns the surface temperature of the earth, computed from the state u.
"""
function get_land_temp_from_state(land_sim, u)
    # required by viz_explorer.jl
    return ClimaLSM.surface_temperature(land_sim.model, u, land_sim.integrator.p, land_sim.integrator.t)
end

"""
    dss_state!(sim::BucketSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on bucket land simulations.
The `dss!` function of ClimaLSM must be called because it uses either the 2D
or 3D dss buffer stored in the cache depending on space of each variable in
`sim.integrator.u`.
"""
function dss_state!(sim::BucketSimulation)
    ClimaLSM.dss!(sim.integrator.u, sim.integrator.p, sim.integrator.t)
end
