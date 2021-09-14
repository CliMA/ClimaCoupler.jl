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

function ∑tendencies_atm_!(dY, Y, (parameters, T_sfc), t)
    println("minimal")
    return dY
end


zero_3vector(el) = Geometry.Cartesian3Vector(el)

function ∑tendencies_atm!(dY, Y, (parameters, T_sfc), t)

    UnPack.@unpack Ch, Cd, f, ν, ug, vg, C_p, MSLP, R_d, R_m, C_v, grav = parameters

    @show Y.x
    (Yc, Yf, F_sfc) = Y.x
    (dYc, dYf, dF_sfc) = dY.x

    UnPack.@unpack ρ, u, v, ρθ = Yc
    UnPack.@unpack w = Yf
    dρ = dYc.ρ
    du = dYc.u
    dv = dYc.v    
    dρθ = dYc.ρθ
    dw = dYf.w
  
    
    #=space_c, space_f = make_function_space(model.domain)
    local_geometry_c = Fields.local_geometry_field(space_c)
    local_geometry_f = Fields.local_geometry_field(space_f)

    # functions that make zeros for this model
    zero_scalar(lg) = zero(FT)
    zero_12vector(lg) = Geometry.Cartesian12Vector(zero(FT), zero(FT))
    zero_3vector(lg) = Geometry.Cartesian3Vector(zero(FT))

    ρ = zero_scalar.(local_geometry_c)
    uv = zero_12vector.(local_geometry_c)
    w = zero_3vector.(local_geometry_f) # faces
    ρθ = zero_scalar.(local_geometry_c)=#

    # u = Operators.getidx(uv, Operators.Interior(), 1)
    # v = Operators.getidx(uv, Operators.Interior(), 2)

    # du = Operators.getidx(duv, Operators.Interior(), 1)
    # dv = Operators.getidx(duv, Operators.Interior(), 2)

    # Auxiliary calculations
    u_1 = parent(u)[1]
    v_1 = parent(v)[1]
    ρ_1 = parent(ρ)[1]
    ρθ_1 = parent(ρθ)[1]
    u_wind = sqrt(u_1^2 + v_1^2)

    # surface flux calculations 
    #surface_flux_ρθ = - calculate_sfc_fluxes_energy(DryBulkFormulaWithRadiation(), parameters, T_sfc[1], parent(ρθ)[1] / parent(ρ)[1] , u_1, v_1, ρ_1, t ) ./ C_p
    
    surface_flux_ρθ = - calculate_sfc_fluxes_energy(DryMonin(), parameters, T_sfc[1], parent(ρθ)[1] / parent(ρ)[1] , u_1, v_1, ρ_1, t ) ./ C_p
    surface_flux_u =  - Cd * u_1 * sqrt(u_1^2 + v_1^2)
    surface_flux_v =  - Cd * v_1 * sqrt(u_1^2 + v_1^2)

    # accumulate in the required right units
    @inbounds begin
        dY.x[3][1] = - ρ_1 * surface_flux_u  # 
        dY.x[3][2] = - ρ_1 * surface_flux_v  # 
        dY.x[3][3] = - C_p * surface_flux_ρθ # W / m^2
    end

    # @inbounds begin
    #     dY.x[3][1] = - 10.0 
    #     dY.x[3][2] = - 1.0  
    #     dY.x[3][3] = - 1.0 
    # end

    # Density tendency (located at cell centers)
    # density
    If = Operators.InterpolateC2F()
    ∂f = Operators.GradientC2F()
    ∂c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT))),
        top = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT))),
    )

    @. dρ = -∂c(zero_3vector.(w) * If(ρ))

    # potential temperature
    bottom_flux_ρθ = zero(FT)
    @show surface_flux_ρθ
    bottom_flux_ρθ = surface_flux_ρθ

    If = Operators.InterpolateC2F()
    ∂f = Operators.GradientC2F()
    ∂c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(bottom_flux_ρθ)),
        top = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT))),
    )
    # TODO!: Undesirable casting to vector required
    @. dρθ =
        -∂c(zero_3vector.(w) * If(ρθ)) + ρ * ∂c(Geometry.CartesianVector(ν * ∂f(ρθ / ρ)))

    # # u velocity tendency (located at cell centers)
    # gradc2f = Operators.GradientC2F(top = Operators.SetValue(ug)) # Eq. 4.18
    # gradf2c = Operators.GradientF2C(bottom = Operators.SetValue(Cd * u_wind * u_1)) # Eq. 4.16
    
    # A = Operators.AdvectionC2C(bottom = Operators.SetValue(0.0), top = Operators.SetValue(0.0))
    # @. du = gradf2c(ν * gradc2f(u)) + f * (v - vg) - A(w, u) # Eq. 4.8

    # # v velocity (centers)
    # gradc2f = Operators.GradientC2F(top = Operators.SetValue(vg)) # Eq. 4.18
    # gradf2c = Operators.GradientF2C(bottom = Operators.SetValue(Cd * u_wind * v_1)) # Eq. 4.16

    # A = Operators.AdvectionC2C(bottom = Operators.SetValue(0.0), top = Operators.SetValue(0.0))
    # @. dv = gradf2c(ν * gradc2f(v)) - f * (u - ug) - A(w, v) # Eq. 4.9

    # # w velocity (faces)
    # gradc2f = Operators.GradientC2F()
    # gradf2c = Operators.GradientF2C(bottom = Operators.SetValue(0.0), top = Operators.SetValue(0.0))

    # B = Operators.SetBoundaryOperator(bottom = Operators.SetValue(0.0), top = Operators.SetValue(0.0))
    # If = Operators.InterpolateC2F(bottom = Operators.Extrapolate(), top = Operators.Extrapolate())
    
    # Π(ρθ) = C_p .* (R_d .* ρθ ./ MSLP).^(R_m ./ C_v)
    
    # @. dw = B( -(If(ρθ / ρ) * gradc2f(Π(ρθ))) - grav + gradc2f(ν * gradf2c(w)) - w * If(gradf2c(w))) # Eq. 4.10 # this makes everything unstable... use new ClimaAtmos rhs!
    
    return dY
end


""" Initialize fields located at cell centers in the vertical. """
# function init_centers(zc, parameters)
#     UnPack.@unpack grav, C_p, MSLP, R_d = parameters

#     # temperature
#     Γ = grav / C_p
#     T = max(land_surface_temperature - Γ * zc, minimum_reference_temperature)

#     # pressure
#     p = MSLP * (T / land_surface_temperature)^(grav / (R_d * Γ))

#     if T == minimum_reference_temperature
#         z_top = (land_surface_temperature - minimum_reference_temperature) / Γ
#         H_min = R_d * minimum_reference_temperature / grav
#         p *= exp(-(zc - z_top) / H_min)
#     end

#     # potential temperature
#     θ = land_surface_temperature

#     # density
#     ρ = p / (R_d * θ * (p / MSLP)^(R_d / C_p))

#     # velocties
#     u = 1.0
#     v = 0.0
    
#     return (ρ = ρ, u = u, v = v, ρθ = ρ * θ)
# end
# function init_centers_and_faces(zc)
#     p = parameters
#     grav, C_p, MSLP, R_d, T_surf, T_min_ref, u0, v0, w0 = p.grav, p.C_p, p.MSLP, p.R_d, p.T_surf, p.T_min_ref, p.u0, p.v0, p.w0

#     # temperature
#     Γ = grav / C_p
#     T = max.(land_surface_temperature .- Γ.* zc, minimum_reference_temperature)

#     # pressure
#     p = MSLP .* (T ./ land_surface_temperature) .^(grav ./ (R_d .* Γ))

#     if T == minimum_reference_temperature
#         z_top = (land_surface_temperature .- minimum_reference_temperature) ./ Γ
#         H_min = R_d .* minimum_reference_temperature ./ grav
#         p *= exp.(-(zc .- z_top) ./ H_min)
#     end

#     # potential temperature
#     θ = land_surface_temperature

#     # density
#     ρ = p ./ (R_d * θ .* (p ./ MSLP) .^(R_d ./ C_p))

#     # velocties
#     u = 1.0 .+ zc .* 0.0 
#     v = 0.0 .+ zc .* 0.0 
    
#     return (; ρ = ρ .+ zc .* 0.0 , u = u, v = v , ρθ = ρ .* θ .+ zc .* 0.0,  w = 0.0 .* zf ) 
# end



""" Initialize fields located at cell interfaces in the vertical. """
function init_faces(zf, parameters)
    return (; w = 0.0 .* zf)
end

function init_centers(zc, parameters )
    p = parameters # import properly
    grav, C_p, MSLP, R_d, T_surf, T_min_ref, u0, v0, w0 = p.grav, p.C_p, p.MSLP, p.R_d, p.T_surf, p.T_min_ref, p.u0, p.v0, p.w0

    minimum_reference_temperature =parameters.T_min_ref
    # temperature
    Γ = grav / C_p
    T = max.(T_surf .+ 10.0 .- Γ.* zc, minimum_reference_temperature)

    # pressure
    p = MSLP .* (T ./ T_surf) .^(grav ./ (R_d .* Γ))

    if T == minimum_reference_temperature
        z_top = (T_surf .- minimum_reference_temperature) ./ Γ
        H_min = R_d .* minimum_reference_temperature ./ grav
        p *= exp.(-(zc .- z_top) ./ H_min)
    end

    # potential temperature
    θ = T_surf + 10.0

    # density
    ρ = p ./ (R_d .* θ .* (p ./ MSLP) .^(R_d ./ C_p))
    @show ρ
    @show θ

    # velocties
    u = 1.0 .+ zc .* 0.0 
    v = 0.0 .+ zc .* 0.0 
    @show u
    return (; ρ = ρ .+ zc .* 0.0 , u = u, v = v , ρθ = ρ .* θ .+ zc .* 0.0 ) 
end