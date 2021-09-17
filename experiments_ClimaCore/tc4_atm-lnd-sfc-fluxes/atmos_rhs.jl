# Atmos RHS

"""
This is a slightly modified ClimaAtmos Ekman column model
"""

#=
Ekman column:
    ∂_t ρ =  ∇ (μ ∇ ρ - w ρ) 
    ∂_t ρθ =  ∇ (μ ∇ ρθ - w ρθ)
    ∂_t u =  ∇ (μ ∇ u - w u) + (v - v_g)
    ∂_t v =  ∇ (μ ∇ v - w v) - (u - u_g) 
    ∂_t w =  ∇ (μ ∇ w - w w) - g - c_p θ ∂_z Π  

where 
    Π = (p/p_0)^{R/c_p}

top BCs are insulating and impenetrable:
    ∂_t T = 0           
    u = u_g             
    v = v_g             
    w = 0.0            
    ∂_t ρ = 0     

and bottom BCs use bulk formulae for surface fluxes of heat and momentum:
    ∂_t ρθ = F₃ = -Ch ρ ||u|| (T_sfc - ρθ / ρ)       
    ∂_t u  = F₁ = -Cd u ||u||                        
    ∂_t v  = F₂ = -Cd v ||u||                        
    w = 0.0                                     
    ∂_t ρ = 0                                   

We also use this model to accumulate fluxes it calculates
    ∂_t F_accum = -(F₁, F₂, F₃)
=#

function ∑tendencies_atm!(dY, Y, (parameters, T_sfc), t)
    @unpack Cd, f, ν, uvg, C_p, MSLP, R_d, R_m, C_v, grav = parameters

    # unpack tendencies and state
    (Yc, Yf, F_sfc) = Y.x
    (dYc, dYf, dF_sfc) = dY.x

    UnPack.@unpack ρ, uv, ρθ = Yc
    
    w = Yf
    dρ = dYc.ρ
    duv = dYc.uv
    dρθ = dYc.ρθ
    dw = dYf

    # auxiliary calculations
    ρ_1 = parent(ρ)[1]
    ρθ_1 = parent(ρθ)[1]
    uv_1 = Operators.getidx(uv, Operators.Interior(), 1)
    u_wind = LinearAlgebra.norm(uv_1)

    # surface flux calculations 
    fluxes = calculate_sfc_fluxes(DryMonin(), parameters, T_sfc[1], ρθ_1 / ρ_1, uv_1, ρ_1, t ) 

    surface_flux_ρθ = fluxes.SH
    surface_flux_uv = fluxes.τ
    
    # accumulate in the required right units
    @inbounds begin
        dY.x[3][1] = (surface_flux_ρθ)  
        dY.x[3][2] = (surface_flux_uv)[1]
        dY.x[3][2] = (surface_flux_uv)[2]   
    end
    
    # boundary conditions
    bcs_bottom_uv_flux = Operators.SetValue(surface_flux_uv)
    bcs_top_uv_value = Operators.SetValue(uvg)
 
    bcs_bottom_ρθ_flux = Operators.SetValue(Geometry.Cartesian3Vector(surface_flux_ρθ))
    bcs_top_ρθ_flux = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT)))

    bcs_bottom_ρ_flux = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT)))
    bcs_top_ρ_flux = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT)))

    # density
    If = Operators.InterpolateC2F()
    ∂f = Operators.GradientC2F()
    ∂c = Operators.DivergenceF2C(
        bottom = bcs_bottom_ρ_flux ,
        top = bcs_top_ρ_flux ,
    )
    @. dρ = -∂c(w * If(ρ))

    # potential temperature
    If = Operators.InterpolateC2F()
    ∂f = Operators.GradientC2F()
    ∂c = Operators.DivergenceF2C(
        bottom = bcs_bottom_ρθ_flux,
        top = bcs_top_ρθ_flux,
    )
    # TODO!: Undesirable casting to vector required
    @. dρθ =
        -∂c(w * If(ρθ)) + ρ * ∂c(Geometry.CartesianVector(ν * ∂f(ρθ / ρ)))

    A = Operators.AdvectionC2C(
        bottom = Operators.SetValue(Geometry.Cartesian12Vector(0.0, 0.0)),
        top = Operators.SetValue(Geometry.Cartesian12Vector(0.0, 0.0)),
    )

    # uv
    ∂c = Operators.DivergenceF2C(bottom = bcs_bottom_uv_flux)
    ∂f = Operators.GradientC2F(top = bcs_top_uv_value)
    duv .= (uv .- Ref(uvg)) .× Ref(Geometry.Cartesian3Vector(f))

    @. duv += ∂c(ν * ∂f(uv)) - A(w, uv)

    # w
    If = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    ∂f = Operators.GradientC2F()
    ∂c = Operators.GradientF2C()
    Af = Operators.AdvectionF2F()
    divf = Operators.DivergenceC2F()
    B = Operators.SetBoundaryOperator(
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT))),
        top = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT))),
    )
    Φ(z) = grav * z
    Π(ρθ) = C_p * (R_d * ρθ / MSLP)^(R_m / C_v)
    zc = Fields.coordinate_field(axes(ρ))
    @. dw = B(
        Geometry.CartesianVector(-(If(ρθ / ρ) * ∂f(Π(ρθ))) - ∂f(Φ(zc))) + divf(ν * ∂c(w)) - Af(w, w),
    )

    return dY
end


""" 
Initialize fields located at cell centers in the vertical. 
"""
function init_ekman_column_1d_c(z, params)
    @unpack grav, C_p, MSLP, R_d, T_surf_atm, T_min_ref, u0, v0, w0 = params

    T_surf = T_surf_atm
    Γ = grav / C_p
    T = max(T_surf - Γ * z, T_min_ref)
    p = MSLP * (T / T_surf)^(grav / (R_d * Γ))
    if T == T_min_ref
        z_top = (T_surf - T_min_ref) / Γ
        H_min = R_d * T_min_ref / grav
        p *= exp(-(z - z_top) / H_min)
    end
    θ = T_surf # potential temperature

    ρ = p / (R_d * θ * (p / MSLP)^(R_d / C_p))

    # velocity
    uv= Geometry.Cartesian12Vector(u0, v0) # u, v components

    # potential temperature
    ρθ = ρ * T_surf

    return (ρ = ρ, uv = uv, ρθ = ρθ)
end


""" 
Initialize fields located at cell interfaces in the vertical. 
"""
function init_ekman_column_1d_f(z, params)
    @unpack grav, C_p, MSLP, R_d, T_surf_atm, T_min_ref, u0, v0, w0 = params

    w = Geometry.Cartesian3Vector(w0) # w component

    return (w = w)
end
