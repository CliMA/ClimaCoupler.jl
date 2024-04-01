"""
    SurfaceStub

On object containing simulation-like info, used as a stub or for prescribed data.
"""
struct SurfaceStub{I} <: SurfaceModelSimulation
    cache::I
end

"""
    stub_init(cache)

Initialization function for SurfaceStub simulation type.
"""
stub_init(cache) = SurfaceStub(cache)

## Extensions of Interfacer.jl functions

"""
    get_field(::SurfaceStub, ::Val)

A getter function, that should not allocate. If undefined, it returns a descriptive error.
"""
get_field(sim::SurfaceStub, ::Val{:air_density}) = sim.cache.ρ_sfc
get_field(sim::SurfaceStub, ::Val{:area_fraction}) = sim.cache.area_fraction
get_field(sim::SurfaceStub, ::Val{:beta}) = sim.cache.beta
get_field(sim::SurfaceStub, ::Val{:energy}) = nothing
get_field(sim::SurfaceStub, ::Val{:roughness_buoyancy}) = sim.cache.z0b
get_field(sim::SurfaceStub, ::Val{:roughness_momentum}) = sim.cache.z0m
get_field(sim::SurfaceStub, ::Val{:surface_direct_albedo}) = sim.cache.α_direct
get_field(sim::SurfaceStub, ::Val{:surface_diffuse_albedo}) = sim.cache.α_diffuse
get_field(sim::SurfaceStub, ::Val{:surface_humidity}) =
    TD.q_vap_saturation_generic.(sim.cache.thermo_params, sim.cache.T_sfc, sim.cache.ρ_sfc, sim.cache.phase)
get_field(sim::SurfaceStub, ::Val{:surface_temperature}) = sim.cache.T_sfc
get_field(sim::SurfaceStub, ::Val{:water}) = nothing

"""
    update_field!(sim::SurfaceStub, ::Val{:area_fraction}, field::Fields.Field)

Updates the specified value in the cache of `SurfaceStub`.
"""
function update_field!(sim::SurfaceStub, ::Val{:area_fraction}, field::Fields.Field)
    sim.cache.area_fraction .= field
end
function update_field!(sim::SurfaceStub, ::Val{:surface_temperature}, field::Fields.Field)
    sim.cache.T_sfc .= field
end
function update_field!(sim::SurfaceStub, ::Val{:air_density}, field)
    parent(sim.cache.ρ_sfc) .= parent(field)
end
function update_field!(sim::SurfaceStub, ::Val{:surface_direct_albedo}, field::Fields.Field)
    sim.cache.α_direct .= field
end
function update_field!(sim::SurfaceStub, ::Val{:surface_diffuse_albedo}, field::Fields.Field)
    sim.cache.α_diffuse .= field
end

"""
    name(::ComponentModelSimulation)

Returns simulation name, if defined, or `Unnamed` if not.
"""
name(::SurfaceStub) = "SurfaceStub"


## Extensions of FieldExchanger.jl functions

"""
    update_sim!(::SurfaceStub, csf, area_fraction)

The stub surface simulation only updates the air density (needed for the turbulent flux calculation).
"""
function update_sim!(sim::SurfaceStub, csf, area_fraction)
    update_field!(sim, Val(:air_density), csf.ρ_sfc)
end

"""
    reinit!(cs::SurfaceStub)

The stub surface simulation is not updated by this function. Extends `SciMLBase.reinit!`.
"""
reinit!(::SurfaceStub) = nothing

"""
    step!(::SurfaceStub, t)

The stub surface simulation is not updated by this function. Extends `SciMLBase.step!`.
"""
step!(::SurfaceStub, _) = nothing


## Extensions of FluxCalculator.jl functions

update_turbulent_fluxes_point!(sim::SurfaceStub, fields::NamedTuple, colidx::Fields.ColumnIndex) = nothing
