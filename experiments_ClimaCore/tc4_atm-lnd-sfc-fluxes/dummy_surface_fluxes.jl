export calculate_sfc_fluxes_energy
export DryBulkFormulaWithRadiation, DryBulkFormula, LinearRelaxation, DryMonin

"""
Options for calulating surface fluxes
"""
abstract type SurfaceFluxType end

struct LinearRelaxation <: SurfaceFluxType end
struct DryBulkFormula <: SurfaceFluxType end
struct DryBulkFormulaWithRadiation <: SurfaceFluxType end
struct DryMonin <: SurfaceFluxType end

calculate_sfc_fluxes_energy(formulation::LinearRelaxation, p, T_sfc, T1) = p.λ .* (T_sfc .- T1)

"""
    calculate_sfc_fluxes_energy(formulation::DryBulkFormula, parameters, T_sfc, θ_1, u_1, v_1, ρ_1, t )

- calculates momentum and sensible heat fluxes using a bulk formula with constant drag/transfer coefficients
"""
calculate_sfc_fluxes_energy(formulation::DryBulkFormula, p, T_sfc, T1, u_1, v_1, ρ_1 ) = p.Ch * p.C_p * ρ_1 * sqrt(u_1^2 + v_1^2) * (T_sfc .- T1)


"""
    calculate_sfc_fluxes_energy(formulation::DryBulkFormulaWithRadiation, parameters, T_sfc, θ_1, u_1, v_1, ρ_1, t )

- calculates momentum and sensible heat fluxes using a bulk formula with constant drag/transfer coefficients
- includes a simple representation of surface radiative fluxes with a temporal dependency
"""
function calculate_sfc_fluxes_energy(formulation::DryBulkFormulaWithRadiation, parameters, T_sfc, θ_1, u_1, v_1, ρ_1, t )

    R_LW_incoming = 0.0
    p = parameters
    R_SW = (1-p.α) * p.τ * p.F_sol * (1 .+ sin(t * 2π / p.τ_d) )
    R_LW = p.ϵ * (p.σ * T_sfc .^ 4 ) - R_LW_incoming
    SH   = p.Ch * p.C_p * ρ_1 * sqrt(u_1^2 + v_1^2) * (T_sfc - θ_1) # convert to theta!!
    #p_sfc = p.pₒ
    #LH   = p.lambda * p.g_w * (q_sat(T_sfc, p_sfc) - q_a)

    F_tot = - (R_SW - R_LW - SH )#- LH
end


"""
    calculate_sfc_fluxes_energy(formulation::DryMonin, parameters, T_sfc, θ_1, windspeed_1, ρ_1, t )

- calculates momentum and sensible heat fluxes using drag/transfer coefficients calculated in CliMA's `SurfaceFluxes.jl` module. This uses the Monin Obukhov theory and Nishizawa and Kitamura (2018) method for the finite volume discretization
"""
function calculate_sfc_fluxes(formulation::DryMonin, parameters, T_sfc, θ_1, uv_1, ρ_1, t )

    domain_atm  = Domains.IntervalDomain(0, parameters.atmos_Lz, x3boundary = (:bottom, :top))
    mesh_atm = Meshes.IntervalMesh(domain_atm, nelems = parameters.atmos_Nz) 
    center_space_atm = Spaces.CenterFiniteDifferenceSpace(mesh_atm) 
    z_centers = Fields.coordinate_field(center_space_atm)
    z_in = parent(z_centers)[1]
    
    windspeed_1 = LinearAlgebra.norm(uv_1)
    p = parameters

    z_0m = 1e-5 # eventually DryMonin.z_0m? (use cases: where atmos will get it with prescribed land, or vice versa, or coupled)
    z_0c = z_0m # z_0c = function(z0m, ustar...) eventually
    z_0 = [z_0m, z_0c]
    x_in = [windspeed_1, θ_1] 
    
    x_s = [windspeed_1* 0.0 , T_sfc] # we are assuming T_sfc = T(-\Delta z).
    θ_basic = deepcopy(θ_1) # for buoyancy calculation

    # Initial guesses for MO parameters, these should be a function of state.
    LMO_init = z_in * 100 # Initial value so that ξ_init<<1
    u_star_init = 0.1 * windspeed_1
    th_star_init = 1.0 * θ_1
    MO_param_guess = [LMO_init, u_star_init, th_star_init] # guess based on current atmos state

    output = SF.surface_conditions(
        CLIMAparam_set,
        MO_param_guess,
        x_in,
        x_s,
        z_0,
        θ_basic,
        z_in,
        SF.FVScheme(),
        maxiter = 500, # often unconverged + coeffs v small, check with FMS (esp Mo_params_guess), update for conditional states
    )

    C_exchange = output.C_exchange
    Cd = C_exchange[1] # Cd
    Ch = C_exchange[2] # Ch
    
    # upward sensible potential temperature flux
    SH   = Ch * ρ_1 * windspeed_1 * (T_sfc - θ_1) 

    # wind stress
    τ = Geometry.Cartesian3Vector(Cd * windspeed_1) ⊗ uv_1

    fluxes = (τ = τ, SH = SH, )
end


