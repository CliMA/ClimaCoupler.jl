# Surface flux calculations

"""
    calculate_land_sfc_fluxes(model::ModelSetup, state, aux, t)
- calculate furface fluxes using the bulk gradient diffusion theory
"""
function calculate_land_sfc_fluxes(model::ModelSetup, state, aux, t; MO_params = nothing) # should pass in the coupler state (also move to coupler), so can access states of both models derectly -e.g. callback?
    
    FT = eltype(state)
    p = model.parameters

    p_sfc = p.pₒ
    θ_a = state.ρθ / state.ρ 
    T_a = θ_a
    T_sfc = aux.T_sfc 
    θ_sfc = T_sfc  

    R_SW = (FT(1)-p.α) * p.τ * p.F_sol * (FT(1) .+ sin(t * 2π / p.τ_d) )
    R_LW = p.ϵ * (p.σ * θ_sfc .^ 4 - p.F_a)
    SH   = state.ρ * p.cp_d * p.g_a * (T_sfc - T_a) # g_a could be substituded by g_a=f(MO_params) (see Bonan P90), but first need to modify SurfaceFluxes.jl
      
    #LH   = p.lambda * p.g_w * (q_sat(T_sfc, p_sfc) - q_a) 
    F_tot = - (R_SW - R_LW - SH )#- LH 

end