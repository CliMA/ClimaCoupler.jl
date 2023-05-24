# atmos_init: for ClimaAtmos pre-AMIP interface
import ClimaAtmos
import StaticArrays

driver_file = joinpath(pkgdir(ClimaAtmos), "examples", "hybrid", "driver.jl")
ENV["CI_PERF_SKIP_RUN"] = true
try
    include(driver_file)
catch err
    if err.error !== :exit_profile
        rethrow(err.error)
    end
end
# the clima atmos `integrator` is now defined
struct ClimaAtmosSimulation{F, P, Y, D, I} <: AtmosModelSimulation
    FT::F
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function atmos_init(::Type{FT}, Y, integrator; params = nothing) where {FT}
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.w)
    spaces = (; center_space = center_space, face_space = face_space)
    if :ρe_int in propertynames(Y.c)
        @warn("Running with ρe_int in coupled mode is not tested yet.")
    end

    isnothing(integrator.p.radiation_model) ? error("ClimaAtmos: radiative model not defined") : nothing
    @assert(isnothing(integrator.p.surface_scheme) == false , "No surface scheme selected")
    @assert(integrator.p.atmos.vert_diff isa CA.VerticalDiffusion, "Atmos fluxes not applied via vertical diffusion")

    ClimaAtmosSimulation(FT, params, Y, spaces, integrator)
end

function update_calculated_fluxes_point!(sim::ClimaAtmosSimulation, fields, colidx)
    (; F_ρτxz , F_ρτyz , F_shf , F_lhf , F_evap ) = fields

    ρ_int = Spaces.level(sim.integrator.u.c.ρ , 1)

    @. sim.integrator.p.dif_flux_energy_bc[colidx] = - Geometry.WVector(F_shf[colidx] + F_lhf[colidx]) # Geometry.WVector(outputs[colidx].F_shf + outputs[colidx].F_lhf,)
    @. sim.integrator.p.dif_flux_ρq_tot_bc[colidx] = - Geometry.WVector(F_evap[colidx])
    @. sim.integrator.p.dif_flux_uₕ_bc[colidx] = - Geometry.Contravariant3Vector(sim.integrator.p.surface_normal[colidx]) ⊗ Geometry.Covariant12Vector(Geometry.UVVector(F_ρτxz / ρ_int, F_ρτyz / ρ_int)[colidx],)

end

function update!(sim::ClimaAtmosSimulation, ::Val{:F_evapnergy}, field)
    @. sim.integrator.p.dif_flux_energy_bc = - Geometry.WVector(field)
end

function update!(sim::ClimaAtmosSimulation, ::Val{:F_evapvaporation}, field)
    @. sim.integrator.p.dif_flux_ρq_tot_bc  = - Geometry.WVector(field)
end

function update!(sim::ClimaAtmosSimulation, ::Val{:F_drag}, (F_ρτxz, F_ρτyz))
    ρ_int = Spaces.level(sim.integrator.u.c.ρ , 1)
    surface_normal = Geometry.WVector.(ones(axes(Fields.level(sim.integrator.u.c, 1))))

    @. sim.integrator.p.dif_flux_uₕ_bc  = - Geometry.Contravariant3Vector(sim.integrator.p.surface_normal[colidx]) ⊗ Geometry.Covariant12Vector(Geometry.UVVector(F_ρτxz / ρ_int, F_ρτyz / ρ_int),)
end

function update!(sim::ClimaAtmosSimulation, ::Val{:T_sfc}, field)
    sim.integrator.p.radiation_model.surface_temperature .= RRTMGPI.field2array(field)
end

function update!(sim::ClimaAtmosSimulation, ::Val{:albedo}, field)

    sim.integrator.p.radiation_model.diffuse_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
    sim.integrator.p.radiation_model.direct_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
end


get_net_surface_radiation(sim::ClimaAtmosSimulation) = level(sim.integrator.p.ᶠradiation_flux, half)
get_liquid_precipitation(sim::ClimaAtmosSimulation) = sim.integrator.p.col_integrated_rain
get_snow_precipitation(sim::ClimaAtmosSimulation) = sim.integrator.p.col_integrated_snow

get_height_int_point(sim::ClimaAtmosSimulation, colidx) = Spaces.level(Fields.coordinate_field(sim.integrator.u.c).z, 1)[colidx]
get_height_sfc_point(sim::ClimaAtmosSimulation, colidx) = Spaces.level(Fields.coordinate_field(sim.integrator.u.f).z, half)[colidx]

function get_uv_int_point(sim::ClimaAtmosSimulation, colidx)
    uₕ_int = Geometry.UVVector.(Spaces.level(sim.integrator.u.c.uₕ, 1))[colidx]
    return @. StaticArrays.SVector(uₕ_int.components.data.:1, uₕ_int.components.data.:2)
end

get_thermo_state_point(sim::ClimaAtmosSimulation, colidx)  = Spaces.level(sim.integrator.p.ᶜts[colidx], 1)
get_thermo_params(sim::ClimaAtmosSimulation) = CAP.thermodynamics_params(sim.integrator.p.params)

get_air_density(::ClimaAtmosSimulation, thermo_params, thermo_state) = TD.air_density.(thermo_params, thermo_state)
get_air_temperature(::ClimaAtmosSimulation, thermo_params, thermo_state_int) = TD.air_temperature.(thermo_params, thermo_state_int)
get_cv_m(::ClimaAtmosSimulation, thermo_params, thermo_state_int) = TD.cv_m.(thermo_params, thermo_state_int)
get_gas_constant_air(::ClimaAtmosSimulation, thermo_params, thermo_state_int)  = TD.gas_constant_air.(thermo_params, thermo_state_int)
get_q_vap_saturation_generic(::ClimaAtmosSimulation, thermo_params, thermo_state_int) = TD.q_vap_saturation_generic.(thermo_params, T_sfc, ρ_sfc, TD.Liquid())

get_surface_params(sim::ClimaAtmosSimulation) = CAP.surface_fluxes_params(sim.integrator.p.params)
function get_surface_scheme(sim::ClimaAtmosSimulation)
    if sim.integrator.p.surface_scheme isa CA.MoninObukhovSurface
        return MoninObukhovScheme()
    elseif sim.integrator.p.surface_scheme isa CA.BulkSurface
        return BulkScheme()
    else
        return nothing
    end
end

reinit!(sim::ClimaAtmosSimulation) = reinit!(sim.integrator)
