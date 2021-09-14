export calculate_sfc_fluxes_energy
export DryBulkFormulaWithRadiation, DryBulkFormula, LinearRelaxation, DryMonin

#using Main.SurfaceFluxes: surface_conditions

"""
dummy for SurfaceFluxes.jl
"""
abstract type SurfaceFluxType end

struct LinearRelaxation <: SurfaceFluxType end
struct DryBulkFormula <: SurfaceFluxType end
struct DryBulkFormulaWithRadiation <: SurfaceFluxType end
struct DryMonin <: SurfaceFluxType end

calculate_sfc_fluxes_energy(formulation::LinearRelaxation, p, T_sfc, T1) = p.λ .* (T_sfc .- T1)

calculate_sfc_fluxes_energy(formulation::DryBulkFormula, p, T_sfc, T1, u_1, v_1, ρ_1 ) = p.Ch * p.C_p * ρ_1 * sqrt(u_1^2 + v_1^2) * (T_sfc .- T1)

function calculate_sfc_fluxes_energy(formulation::DryBulkFormulaWithRadiation, parameters, T_sfc, T_1, u_1, v_1, ρ_1, t )

    R_LW_incoming = 0.0
    p = parameters
    R_SW = (1-p.α) * p.τ * p.F_sol * (1 .+ sin(t * 2π / p.τ_d) )
    R_LW = p.ϵ * (p.σ * T_sfc .^ 4 ) - R_LW_incoming
    SH   = p.Ch * p.C_p * ρ_1 * sqrt(u_1^2 + v_1^2) * (T_sfc - T_1) # convert to theta!!
    #p_sfc = p.pₒ
    #LH   = p.lambda * p.g_w * (q_sat(T_sfc, p_sfc) - q_a)

    F_tot = - (R_SW - R_LW - SH )#- LH
end

Π = 1.0  #(p0/ p)^(R / C_p)

#=
v_1 = 270.0
u_1 = 1.0
T_1 = 0.0
ρ_1 = 1.0
T_sfc = 260.0


wθ_flux_star = nothing

scheme = FVScheme()
θ_scale = θ_basic
universal_func = Businger
=#

function calculate_sfc_fluxes_energy(formulation::DryMonin, parameters, T_sfc, T_1, u_1, v_1, ρ_1, t )

    domain_atm  = Domains.IntervalDomain(0, parameters.atmos_Lz, x3boundary = (:bottom, :top)) # struct
    mesh_atm = Meshes.IntervalMesh(domain_atm, nelems = parameters.atmos_Nz) # struct, allocates face boundaries to 5,6
    center_space_atm = Spaces.CenterFiniteDifferenceSpace(mesh_atm) # collection of the above, discretises space into FD and provides coords
    z_centers = Fields.coordinate_field(center_space_atm)

    θ_1 = T_1 * Π
    windspeed_1 = sqrt(u_1^2 + v_1^2)

    p = parameters
#    MO_param_guess = guess based on current atmos state
    z_0m = 0.01 #eventually DryMonin.z_0m? (use cases: where atmos will get it with prescribed land, or vice versa, or coupled)
    z_0c = z_0m#    z_0c = function(z0m, ustar...) eventually?
    z_0 = [z_0m, z_0c]
    x_in = [windspeed_1, θ_1] # input vector of atmos first level
    @show T_1
    x_s = [windspeed_1* 0.0 , T_sfc] # we are assuming T_sfc = T(-\Delta z).

    θ_basic = deepcopy(T_sfc) # for buoyancy calculation

    ## Initial guesses for MO parameters, these should be a function of state.
    LMO_init = 100 # Initial value so that ξ_init<<1
    u_star_init = 0.1 * windspeed_1
    th_star_init = T_sfc
    MO_param_guess = [LMO_init, u_star_init, th_star_init]

    z_in = parent(z_centers)[1]

    output = SF.surface_conditions(
        CLIMAparam_set,
        MO_param_guess,
        x_in,
        x_s,
        z_0,
        θ_basic,
        z_in,
        SF.FVScheme(),
    )

    C_exchange = output.C_exchange
    Ch = C_exchange[2]

    SH   = Ch * p.C_p * ρ_1 * sqrt(u_1^2 + v_1^2) * (T_sfc - θ_1) #(T_sfc - T_1) for TE

    F_tot =  SH
end



# struct SurfaceFluxConditions{FT, VFT}
#     L_MO::FT
#     wθ_flux_star::FT
#     flux::VFT
#     x_star::VFT
#     C_exchange::VFT
# end