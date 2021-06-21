########
# Set up parameters
########
parameters = (
        xmax = 10000,   # horizontal extent of the SingleStack in the x direction [m]
        ymax = 500,     # horizontal extent of the SingleStack in the y direction [m]
        zmax = 10000,   # vertical extent of the SingleStack [m]
        g    = 9.81,    # gravitation acceleration [m /s^2]
        pₒ   = 1.01325e5,                                               # initial surface pressure [Pa]
        cp_d = (8.3144598 / 28.97e-3) / (2/7),                          # specific heat for ideal gas at constant pressure [J / K / kg]
        cv_d = (8.3144598 / 28.97e-3) / (2/7) - 8.3144598 / 28.97e-3,   # specific heat for ideal gas at constant volume [J / K / kg]
        R_d  = 8.3144598 / 28.97e-3,                                    # gas constant for dry air [J / K / kg]
        κ    = 2.0/7.0, # R_d/cp_d
        cₛ   = 340,     # speed of sound [m / s]
        c_s  = 800,     # specific heat for land (soil)  [J / K / kg]
        κ_s  = 0.0,     # soil conductivity [W / m / K] (set to 0 for energy conservation checks)
        h_s  = 100,     # depth of the modelled soil model [m]
        T_h  = 280,     # temperature of soil at depth h [K]
        ρ_s  = 1500,    # density for land (soil) [kg / m^3]
        τ    = 0.9,     # atmospheric transmissivity 
        α    = 0.5,     # surface albedo    
        σ    = 5.67e-8, # Steffan Boltzmann constant [kg / s^3 / K^4]
        g_a  = 0.06,    # aerodynamic conductance for heat transfer [kg / m^2 / s]
        ϵ    = 0.98,    # broadband emissivity / absorptivity of the surface
        F_a  = 0.0,     # downward LW flux from the atmosphere [W / m^2]
        F_sol = 1361,   # incoming solar TOA radiation [W / m^2]
        τ_d   = 10,     # idealized daily cycle period [s]
    )
 
########
# Set up inital conditions
########

# Land (slab) initial condition
T_sfc₀(p, x, y, z) = p.T_h

# Atmos (single stack) initial conditions
θ₀Atmos(p, x, y, z)         = 300.0 + 10.0 * ( z /p.zmax ) 
π_exnerAtmos(p, x, y, z)    = 1.0 - p.g / (p.cp_d * θ₀Atmos(p, x, y, z) ) * z  
ρ₀Atmos(p, x, y, z)         = p.pₒ / (p.R_d * θ₀Atmos(p, x, y, z)) * (π_exnerAtmos(p, x, y, z))^(p.cv_d/p.R_d)
ρθ₀Atmos(p, x, y, z)        = ρ₀Atmos(p, x, y, z) * θ₀Atmos(p, x, y, z) 
ρu₀Atmos(p, x, y, z)        = @SVector [1,0,0]