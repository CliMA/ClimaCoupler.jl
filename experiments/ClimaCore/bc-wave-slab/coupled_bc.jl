coupler_atmos_boundary_flux(_, ..) = nothing

# momentum
struct CouplerVectorFlux{C, F} <: AbstractBoundary
    coefficients::C
    flux::F
end
function get_boundary_flux(model, bc::CouplerVectorFlux, var::Fields.Field, Y, Ya) # TODO: adapt for Field BC (CC#325)
    flux = Geometry.Covariant3Vector(FT(1)) ⊗ Geometry.Covariant12Vector(FT(1), FT(1)) * parent(bc.flux)[1]
end
"""
get_boundary_flux(
    model,
    bc::CouplerVectorFlux,
    ρc::Fields.Field,
    Y,
    Ya,
)

This is an extenison of ClimaAtmos.get_boundary_flux for momentum
"""
function coupler_atmos_boundary_flux(bc::CouplerVectorFlux, atmos_sim, slab_sim, F_a_space, F_s_space)

    Y_atm = atmos_sim.integrator.u
    FT = eltype(Y_atm)

    uv = Y_atm.base.uh

    uh1_cov = Fields.level(uv.components.data.:1, 1) # TODO: change to Field - pending ClimaCore PR #325 
    uh2_cov = Fields.level(uv.components.data.:2, 1) # TODO: change to Field - pending ClimaCore PR #325 

    # physical scale (wind * coeff)
    uv_1 = Fields.level(Geometry.UVVector.(uv), 1)
    u_wind = norm(uv_1) # TODO: this will need to be spatially variable

    # unit vector in W direction on central space
    normal_unit_vector = Fields.level(similar(Geometry.WVector.(Y_atm.thermodynamics.ρe_tot)), 1) # TODO: for topography, needs to be perp to surface
    parent(normal_unit_vector) .= FT(1)
    local_geometry = Fields.local_geometry_field(axes(normal_unit_vector))

    # contravariant scale
    scale_con =
        Geometry.Contravariant3Vector.(
            Geometry.WVector.(normal_unit_vector .* bc.coefficients .* u_wind),
            local_geometry,
        )

    parent(bc.flux) .=
        parent(Geometry.Contravariant3Vector.(scale_con) .⊗ Geometry.Covariant12Vector.(uh1_cov, uh2_cov))
end

# energy
struct CouplerEnergyFlux{C, F} <: AbstractBoundary
    coefficients::C
    flux::F
end

"""
get_boundary_flux(
    model,
    bc::CouplerEnergyFlux,
    ρc::Fields.Field,
    Y,
    Ya,
)

This is an extenison of ClimaAtmos.get_boundary_flux for enthalpy
"""
function get_boundary_flux(model, bc::CouplerEnergyFlux, ρc::ClimaCore.Fields.Field, Y, Ya)
    bulk_flux = bc.flux
    ClimaCore.Geometry.WVector.(parent(bulk_flux)[1]) #flux = Geometry.Contravariant3Vector.(bulk_flux) # need CC extension - waiting for PR (https://github.com/CliMA/ClimaCore.jl/pull/325)
end

"""
coupler_boundary_flux(
    model,
    bc::CouplerEnergyFlux,
    ρc::Fields.Field,
    Y,
    Ya,
)

Vertical fluxes for arbitrary variables with the bulk aerodynamic turbulent formula,  
given variable exchange coefficients. (e.g. for energy, or tracers)
Currently this is a simplified version of the bulk formula for testing. 
"""
function coupler_atmos_boundary_flux(bc::CouplerEnergyFlux, atmos_sim, slab_sim, F_a_space, F_s_space)

    Y_atm = atmos_sim.integrator.u
    FT = eltype(Y_atm)

    T_1 = ClimaCore.Fields.level(calculate_temperature(atmos_sim, FT), 1)

    ρ_1 = ClimaCore.Fields.level(Y_atm.base.ρ, 1)

    # remap slab > atmos
    T_sfc_atmos_grid = ClimaCore.Fields.zeros(F_a_space)
    T_sfc_slab_grid = copy(slab_sim.integrator.u.T_sfc)
    # ClimaCoreTempestRemap.remap!(T_sfc_atmos_grid, T_sfc_slab_grid, R_slab2atm)
    dummmy_remap!(T_sfc_slab_grid, T_sfc_atmos_grid)

    # save in atmos BCs
    c_pd = CLIMAParameters.Planet.cp_d(atmos_sim.model.parameters)
    # parent(bc.flux) .= bc.coefficients .* c_pd .* parent(ρ_1) .* (parent(T_1) .- parent(T_sfc_atmos_grid)) # parent(bc.flux) .= bc.coefficients .* parent(ρ_1 .* u_wind) .* (parent(c_1) .- parent(c_sfc)) # TODO: make neater - same space but different instance error otherwise
    parent(bc.flux) .= FT(10)
    # TODO: make neater - same space but different instance error otherwise
    nothing
end

# code below should be moved to CA / thermodyn 
using LinearAlgebra: norm_sqr
using Thermodynamics

function calculate_temperature(atmos_sim, FT)
    model = atmos_sim.model
    Y = atmos_sim.integrator.u
    p = calculate_pressure(Y, Y, model.thermodynamics, model.moisture, model.moisture, model.parameters, FT)
    R = FT(287)
    T = p ./ R ./ Y.base.ρ
end

@inline function calculate_gravitational_potential(Y, Ya, params, FT)
    g::FT = CLIMAParameters.Planet.grav(params)
    ρ = Y.base.ρ
    z = Fields.coordinate_field(axes(ρ)).z

    return @. g * z
end

@inline function calculate_pressure(Y, Ya, i, ii, iii, params, FT)
    ρ = Y.base.ρ
    uh = Y.base.uh
    w = Y.base.w
    ρe_tot = Y.thermodynamics.ρe_tot

    interp_f2c = Operators.InterpolateF2C()

    z = Fields.coordinate_field(axes(ρ)).z
    K = calculate_kinetic_energy(Y, Y, params, FT)
    Φ = calculate_gravitational_potential(Y, Y, params, FT)

    e_int = @. ρe_tot / ρ - Φ - K
    p = Thermodynamics.air_pressure.(Thermodynamics.PhaseDry.(params, e_int, ρ))

    return p
end

@inline function calculate_kinetic_energy(Y, Ya, params, FT)
    cuₕ = Y.base.uh # Covariant12Vector on centers
    fw = Y.base.w # Covariant3Vector on faces
    If2c = Operators.InterpolateF2C()

    cuvw = Geometry.Covariant123Vector.(cuₕ) .+ Geometry.Covariant123Vector.(If2c.(fw))
    cK = @. norm_sqr(cuvw) / 2
    return cK
end
