########
# Set up parameters
########
parameters = (
    a        = get_planet_parameter(:planet_radius),
    Ω        = get_planet_parameter(:Omega),
    g        = get_planet_parameter(:grav),
    κ        = get_planet_parameter(:kappa_d),
    R_d      = get_planet_parameter(:R_d),
    R_v      = get_planet_parameter(:R_v),
    cv_d     = get_planet_parameter(:cv_d),
    cv_v     = get_planet_parameter(:cv_v),
    cv_l     = get_planet_parameter(:cv_l),
    cv_i     = get_planet_parameter(:cv_i),
    cp_d     = get_planet_parameter(:cp_d), # specific heat for ideal gas at constant pressure [J / K / kg]
    cp_v     = get_planet_parameter(:cp_v), # specific heat for ideal gas at constant volume [J / K / kg]
    cp_l     = get_planet_parameter(:cp_l),
    cp_i     = get_planet_parameter(:cp_i),
    molmass_ratio = get_planet_parameter(:molmass_dryair)/get_planet_parameter(:molmass_water),
    γ        = get_planet_parameter(:cp_d)/get_planet_parameter(:cv_d),
    pₒ       = get_planet_parameter(:MSLP), # initial surface pressure [Pa]
    pₜᵣ      = get_planet_parameter(:press_triple),
    Tₜᵣ      = get_planet_parameter(:T_triple),
    T_0      = get_planet_parameter(:T_0),
    LH_v0    = get_planet_parameter(:LH_v0), # latent heat of vaporisation at 0C [J / kg]
    e_int_v0 = get_planet_parameter(:e_int_v0),
    e_int_i0 = get_planet_parameter(:e_int_i0),
    H        = 30e3,
    k        = 3.0,
    Γ        = 0.005,
    T_E      = 300,
    T_P      = 271.0,
    b        = 2.0,
    z_t      = 15e3,
    λ_c      = π / 9,
    ϕ_c      = 2 * π / 9,
    V_p      = 1.0,
    ϕ_w      = 2*π/9,
    p_w      = 3.4e4,
    q₀       = 0.0, #0.018,
    qₜ       = 0.0, #1e-12,
    ΔT       = 29.0,
    Tₘᵢₙ     = 271.0,
    Δϕ       = 26π/180.0,
    day      = 86400,
    T_ref    = 255,
    τ_precip = 100.0,
    p0       = 1e5,
    Cₑ       = 0.0005, # bulk transfer coefficient for sensible heat
    Cₗ       = 0.0, # bulk transfer coefficient for latent heat
    Mᵥ       = 0.608,
    c_o = 3.93e3,       # specific heat for ocean  [J / K / kg]
    T_h = 280,      # initial temperature of surface ocean layer [K]
    h_o = 100,      # depth of the modelled ocean mixed layer [m]
    F_sol = 1361,   # incoming solar TOA radiation [W / m^2]
    τ = 0.9,        # atmospheric transmissivity 
    α = 0.5,        # surface albedo    
    σ = 5.67e-8,    # Steffan Boltzmann constant [kg / s^3 / K^4]
    #g_a = 0.06,     # aerodynamic conductance for heat transfer [kg / m^2 / s]
    ϵ = 0.98,       # broadband emissivity / absorptivity
    F_a = 0.0,      # LW flux from the atmosphere [W / m^2]
    ρ_o = 1026,     # density for ocean [kg / m^3]
    diurnal_period = 10, # idealized daily cycle period [s]
)

# Mask to pick out the coupled boundary in the MPIStateArrays (here at altitude = 0 m)
epss = sqrt(eps(Float64))
boundary_mask( 𝒫, xc, yc, zc ) = @. abs(( xc^2 + yc^2 + zc^2 )^0.5 - 𝒫.a - 𝒫.H) < epss

########
# Set up inital conditions
########

# 1. Land (ocean) initial condition
T_sfc₀(𝒫, xc, yc, zc) = 270.0 #T_sfc₀(𝒫, xc, yc, zc) = boundary_mask( 𝒫, xc, yc, zc ) * 𝒫.T_h

# 2. Atmos (single stack) initial conditions
# additional initial condition parameters
T₀(𝒫)   = 0.5 * (𝒫.T_E + 𝒫.T_P)
A(𝒫)    = 1.0 / 𝒫.Γ
B(𝒫)    = (T₀(𝒫) - 𝒫.T_P) / T₀(𝒫) / 𝒫.T_P
C(𝒫)    = 0.5 * (𝒫.k + 2) * (𝒫.T_E - 𝒫.T_P) / 𝒫.T_E / 𝒫.T_P
H(𝒫)    = 𝒫.R_d * T₀(𝒫) / 𝒫.g
d_0(𝒫)  = 𝒫.a / 6

# convenience functions that only depend on height
τ_z_1(𝒫,r)   = exp(𝒫.Γ * (r - 𝒫.a) / T₀(𝒫))
τ_z_2(𝒫,r)   = 1 - 2 * ((r - 𝒫.a) / 𝒫.b / H(𝒫))^2
τ_z_3(𝒫,r)   = exp(-((r - 𝒫.a) / 𝒫.b / H(𝒫))^2)
τ_1(𝒫,r)     = 1 / T₀(𝒫) * τ_z_1(𝒫,r) + B(𝒫) * τ_z_2(𝒫,r) * τ_z_3(𝒫,r)
τ_2(𝒫,r)     = C(𝒫) * τ_z_2(𝒫,r) * τ_z_3(𝒫,r)
τ_int_1(𝒫,r) = A(𝒫) * (τ_z_1(𝒫,r) - 1) + B(𝒫) * (r - 𝒫.a) * τ_z_3(𝒫,r)
τ_int_2(𝒫,r) = C(𝒫) * (r - 𝒫.a) * τ_z_3(𝒫,r)
F_z(𝒫,r)     = (1 - 3 * ((r - 𝒫.a) / 𝒫.z_t)^2 + 2 * ((r - 𝒫.a) / 𝒫.z_t)^3) * ((r - 𝒫.a) ≤ 𝒫.z_t)

# convenience functions that only depend on longitude and latitude
d(𝒫,λ,ϕ)     = 𝒫.a * acos(sin(ϕ) * sin(𝒫.ϕ_c) + cos(ϕ) * cos(𝒫.ϕ_c) * cos(λ - 𝒫.λ_c))
c3(𝒫,λ,ϕ)    = cos(π * d(𝒫,λ,ϕ) / 2 / d_0(𝒫))^3
s1(𝒫,λ,ϕ)    = sin(π * d(𝒫,λ,ϕ) / 2 / d_0(𝒫))
cond(𝒫,λ,ϕ)  = (0 < d(𝒫,λ,ϕ) < d_0(𝒫)) * (d(𝒫,λ,ϕ) != 𝒫.a * π)

# base-state thermodynamic variables
I_T(𝒫,ϕ,r)   = (cos(ϕ) * r / 𝒫.a)^𝒫.k - 𝒫.k / (𝒫.k + 2) * (cos(ϕ) * r / 𝒫.a)^(𝒫.k + 2)
Tᵥ(𝒫,ϕ,r)    = (τ_1(𝒫,r) - τ_2(𝒫,r) * I_T(𝒫,ϕ,r))^(-1) * (𝒫.a/r)^2
p(𝒫,ϕ,r)     = 𝒫.pₒ * exp(-𝒫.g / 𝒫.R_d * (τ_int_1(𝒫,r) - τ_int_2(𝒫,r) * I_T(𝒫,ϕ,r)))
q(𝒫,ϕ,r)     = 0.0 #(p(𝒫,ϕ,r) > 𝒫.p_w) ? 𝒫.q₀ * exp(-(ϕ / 𝒫.ϕ_w)^4) * exp(-((p(𝒫,ϕ,r) - 𝒫.pₒ) / 𝒫.p_w)^2) : 𝒫.qₜ

# base-state velocity variables
U(𝒫,ϕ,r)  = 𝒫.g * 𝒫.k / 𝒫.a * τ_int_2(𝒫,r) * Tᵥ(𝒫,ϕ,r) * ((cos(ϕ) * r / 𝒫.a)^(𝒫.k - 1) - (cos(ϕ) * r / 𝒫.a)^(𝒫.k + 1))
u(𝒫,ϕ,r)  = -𝒫.Ω * r * cos(ϕ) + sqrt((𝒫.Ω * r * cos(ϕ))^2 + r * cos(ϕ) * U(𝒫,ϕ,r))
v(𝒫,ϕ,r)  = 0.0
w(𝒫,ϕ,r)  = 0.0

# velocity perturbations
δu(𝒫,λ,ϕ,r)  = -16 * 𝒫.V_p / 3 / sqrt(3) * F_z(𝒫,r) * c3(𝒫,λ,ϕ) * s1(𝒫,λ,ϕ) * (-sin(𝒫.ϕ_c) * cos(ϕ) + cos(𝒫.ϕ_c) * sin(ϕ) * cos(λ - 𝒫.λ_c)) / sin(d(𝒫,λ,ϕ) / 𝒫.a) * cond(𝒫,λ,ϕ)
δv(𝒫,λ,ϕ,r)  = 16 * 𝒫.V_p / 3 / sqrt(3) * F_z(𝒫,r) * c3(𝒫,λ,ϕ) * s1(𝒫,λ,ϕ) * cos(𝒫.ϕ_c) * sin(λ - 𝒫.λ_c) / sin(d(𝒫,λ,ϕ) / 𝒫.a) * cond(𝒫,λ,ϕ)
δw(𝒫,λ,ϕ,r)  = 0.0

# CliMA prognostic variables
# compute the total energy
uˡᵒⁿ(𝒫,λ,ϕ,r)   = u(𝒫,ϕ,r) + δu(𝒫,λ,ϕ,r)
uˡᵃᵗ(𝒫,λ,ϕ,r)   = v(𝒫,ϕ,r) + δv(𝒫,λ,ϕ,r)
uʳᵃᵈ(𝒫,λ,ϕ,r)   = w(𝒫,ϕ,r) + δw(𝒫,λ,ϕ,r)

# cv_m and R_m for moist experiment
cv_m(𝒫,ϕ,r)  = 𝒫.cv_d + (𝒫.cv_v - 𝒫.cv_d) * q(𝒫,ϕ,r)
R_m(𝒫,ϕ,r) = 𝒫.R_d * (1 + (𝒫.molmass_ratio - 1) * q(𝒫,ϕ,r))

T(𝒫,ϕ,r) = Tᵥ(𝒫,ϕ,r) / (1 + 𝒫.Mᵥ * q(𝒫,ϕ,r)) 
e_int(𝒫,λ,ϕ,r)  = cv_m(𝒫,ϕ,r) * (T(𝒫,ϕ,r) - 𝒫.T_0) + q(𝒫,ϕ,r) * 𝒫.e_int_v0
e_kin(𝒫,λ,ϕ,r)  = 0.5 * ( uˡᵒⁿ(𝒫,λ,ϕ,r)^2 + uˡᵃᵗ(𝒫,λ,ϕ,r)^2 + uʳᵃᵈ(𝒫,λ,ϕ,r)^2 )
e_pot(𝒫,λ,ϕ,r)  = 𝒫.g * r

ρ₀(𝒫,λ,ϕ,r)    = p(𝒫,ϕ,r) / R_m(𝒫,ϕ,r) / T(𝒫,ϕ,r)
ρuˡᵒⁿ(𝒫,λ,ϕ,r) = ρ₀(𝒫,λ,ϕ,r) * uˡᵒⁿ(𝒫,λ,ϕ,r)
ρuˡᵃᵗ(𝒫,λ,ϕ,r) = ρ₀(𝒫,λ,ϕ,r) * uˡᵃᵗ(𝒫,λ,ϕ,r)
ρuʳᵃᵈ(𝒫,λ,ϕ,r) = ρ₀(𝒫,λ,ϕ,r) * uʳᵃᵈ(𝒫,λ,ϕ,r)
ρe(𝒫,λ,ϕ,r) = ρ₀(𝒫,λ,ϕ,r) * (e_int(𝒫,λ,ϕ,r) + e_kin(𝒫,λ,ϕ,r) + e_pot(𝒫,λ,ϕ,r))
ρq(𝒫,λ,ϕ,r) = ρ₀(𝒫,λ,ϕ,r) * q(𝒫,ϕ,r)

# Cartesian Representation (boiler plate really)
ρ₀ᶜᵃʳᵗ(𝒫, x...)  = ρ₀(𝒫, lon(x...), lat(x...), rad(x...))
ρu⃗₀ᶜᵃʳᵗ(𝒫, x...) = (   ρuʳᵃᵈ(𝒫, lon(x...), lat(x...), rad(x...)) * r̂(x...)
                     + ρuˡᵃᵗ(𝒫, lon(x...), lat(x...), rad(x...)) * ϕ̂(x...)
                     + ρuˡᵒⁿ(𝒫, lon(x...), lat(x...), rad(x...)) * λ̂(x...) )
ρeᶜᵃʳᵗ(𝒫, x...) = ρe(𝒫, lon(x...), lat(x...), rad(x...))
ρqᶜᵃʳᵗ(𝒫, x...) = ρq(𝒫, lon(x...), lat(x...), rad(x...))


#####
# Held-Suarez Forcing
#####
struct HeldSuarezForcing{S} <: AbstractPhysicsComponent
    parameters::S
end

FT = Float64
day = 86400
held_suarez_parameters = (;
    k_a = FT(1 / (40 * day)),
    k_f = FT(1 / day),
    k_s = FT(1 / (4 * day)),
    ΔT_y = FT(65),
    Δθ_z = FT(10),
    T_equator = FT(294),
    T_min = FT(200),
    σ_b = FT(7 / 10),
    R_d  = parameters.R_d,
    day  = parameters.day,
    grav = parameters.g,
    cp_d = parameters.cp_d,
    cv_d = parameters.cv_d,
    MSLP = parameters.p0,  
)

######
# Modified Held-Suarez Forcing
######
function calc_component!(
    source,
    hsf::HeldSuarezForcing,
    state,
    aux,
    physics,
)
    FT = eltype(state)
    
    _R_d  = hsf.parameters.R_d
    _day  = hsf.parameters.day
    _grav = hsf.parameters.grav
    _cp_d = hsf.parameters.cp_d
    _cv_d = hsf.parameters.cv_d
    _p0   = hsf.parameters.MSLP  

    # Parameters
    T_ref = FT(255)

    # Extract the state
    ρ = state.ρ
    ρu = state.ρu
    ρe = state.ρe
    Φ = aux.Φ
    
    x = aux.x
    y = aux.y
    z = aux.z
    coord = @SVector[x,y,z]

    p = calc_pressure(physics.eos, state, aux, physics.parameters)
    T = p / (ρ * _R_d)

    # Held-Suarez parameters
    k_a  = hsf.parameters.k_a
    k_f  = hsf.parameters.k_f
    k_s  = hsf.parameters.k_s
    ΔT_y = hsf.parameters.ΔT_y
    Δθ_z = hsf.parameters.Δθ_z
    T_equator = hsf.parameters.T_equator
    T_min = hsf.parameters.T_min
    σ_b = hsf.parameters.σ_b

    # Held-Suarez forcing
    φ = @inbounds asin(coord[3] / norm(coord, 2))

    #TODO: replace _p0 with dynamic surfce pressure in Δσ calculations to account
    #for topography, but leave unchanged for calculations of σ involved in T_equil
    σ = p / _p0
    exner_p = σ^(_R_d / _cp_d)
    Δσ = (σ - σ_b) / (1 - σ_b)
    height_factor = max(0, Δσ)
    T_equil = (T_equator - ΔT_y * sin(φ)^2 - Δθ_z * log(σ) * cos(φ)^2) * exner_p
    T_equil = max(T_min, T_equil)
    k_T = k_a + (k_s - k_a) * height_factor * cos(φ)^4
    k_v = k_f * height_factor

    # horizontal projection
    k = coord / norm(coord)
    P = I - k * k'

    # Apply Held-Suarez forcing
    source.ρu -= k_v * P * ρu
    source.ρe -= k_T * ρ * _cv_d * (T - T_equil)
    return nothing
end



struct FluxAccumulator <: AbstractPhysicsComponent end

function calc_component!(
    source,
    hsf::FluxAccumulator,
    state,
    aux,
    physics,)
    
    E, H = calc_ocean_sfc_fluxes(physics, state, aux) 
    source.F_ρe_accum = - (E + H) # latent + sensible heat fluxes [W/m^2]
end


"""
    calculate_land_sfc_fluxes(model::DryAtmosModel, state, aux)
- calculate furface fluxes using the bulk gradient diffusion theory
"""
function calc_ocean_sfc_fluxes(physics, state⁻, aux⁻; MO_params = nothing) # should pass in the coupler state (also move to coupler), so can access states of both models derectly -e.g. callback?

    # Apply bulks laws using the tangential velocity as energy flux
    ρ = state⁻.ρ
    ρu = state⁻.ρu
    ρq = state⁻.ρq
    eos = physics.eos
    parameters = physics.parameters

    # vertical unit vector
    n̂ = aux⁻.∇Φ / parameters.g

    # obtain surface fields from bcs
    Cₕ = parameters.Cₑ
    Cₗ = parameters.Cₗ
    LH_v0 = parameters.LH_v0
    T_sfc = Float64(280)#aux⁻.T_sfc

    # magnitude of tangential velocity (usually called speed)
    u = ρu / ρ
    speed_tangential =Float64(1)# norm((I - n̂ ⊗ n̂) * u)
        
    # sensible heat flux
    cp = calc_heat_capacity_at_constant_pressure(eos, state⁻, parameters)
    T = calc_air_temperature(eos, state⁻, aux⁻, parameters)
    H = -  Cₕ * speed_tangential * cp * Float64(10)#(T - T_sfc)

    #

    # latent heat flux
    q = ρq / ρ
    q_tot_sfc  = calc_saturation_specific_humidity(ρ, T_sfc, parameters) 
    E =  Float64(0.0) #ρ * Cₗ * speed_tangential * LH_v0 * (q - q_tot_sfc) 

    return E, H # E = 0, H = constant

end