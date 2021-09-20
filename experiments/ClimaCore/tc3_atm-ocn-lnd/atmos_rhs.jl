# Atmos RHS

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

    @unpack Cd, f, ν, uvg, C_p, MSLP, R_d, R_m, C_v, grav = model.parameters

    # unpack tendencies and state
    (Yc, Yf, F_sfc) = Y.x
    (dYc, dYf, dF_sfc) = dY.x

    UnPack.@unpack ρ, uv, ρθ = Yc
    
    w = Yf
    dρ = dYc.ρ
    duv = dYc.uv
    dρθ = dYc.ρθ
    dw = dYf

    # Auxiliary calculations
    ρ_1 = parent(ρ)[1]
    ρθ_1 = parent(ρθ)[1]
    uv_1 = Operators.getidx(uv, Operators.Interior(), 1)
    u_wind = LinearAlgebra.norm(uv_1)

    # surface flux calculations 
    surface_flux_ρθ = - calculate_sfc_fluxes_energy(DryBulkFormulaWithRadiation(), parameters, T_sfc[1], parent(ρθ)[1] / parent(ρ)[1] , u_1, v_1, ρ_1, t ) ./ C_p
    surface_flux_u =  - Cd * u_1 * sqrt(u_1^2 + v_1^2)
    surface_flux_v =  - Cd * v_1 * sqrt(u_1^2 + v_1^2)

    # accumulate in the required right units
    @inbounds begin
        dF_sfc[1] = - ρ_1 * surface_flux_u  # 
        dF_sfc[2] = - ρ_1 * surface_flux_v  # 
        dF_sfc[3] = - C_p * surface_flux_ρθ # W / m^2
    end

    # Density tendency (located at cell centers)
    If = Operators.InterpolateC2F()
    ∂f = Operators.GradientC2F()
    ∂c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
        top = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
    )
    @. dρ = -∂c(w * If(ρ))

    # Potential temperature tendency (located at cell centers)
    # ∂_t ρθ =  ∇ (μ ∇ ρθ - w ρθ)
    If = Operators.InterpolateC2F()
    ∂f = Operators.GradientC2F()
    ∂c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(surface_flux_ρθ)),
        top = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
    )
    # TODO!: Undesirable casting to vector required
    @. dρθ =
        -∂c(w * If(ρθ)) + ρ * ∂c(Geometry.CartesianVector(ν * ∂f(ρθ / ρ)))
 
    # uv velocity tendency (located at cell centers)
    uv_1 = Operators.getidx(uv, Operators.Interior(), 1)
    u_wind = LinearAlgebra.norm(uv_1)

    A = Operators.AdvectionC2C(
        bottom = Operators.SetValue(Geometry.Cartesian12Vector(0.0, 0.0)),
        top = Operators.SetValue(Geometry.Cartesian12Vector(0.0, 0.0)),
    )

    # uv
    bcs_bottom = Operators.SetValue(
        Geometry.Cartesian3Vector(Cd * u_wind) ⊗ uv_1,
    )
    bcs_top = Operators.SetValue(uvg)
    ∂c = Operators.DivergenceF2C(bottom = bcs_bottom)
    ∂f = Operators.GradientC2F(top = bcs_top)
    duv .= (uv .- Ref(uvg)) .× Ref(Geometry.Cartesian3Vector(f))
    @. duv += ∂c(ν * ∂f(uv)) - A(w, uv)
    
    # w velocity (faces)
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