########
# Set up parameters
########
parameters = (
        R_d  = 8.3144598 / 28.97e-3,
        pₒ   = 1.01325e5,
        κ    = 2.0/7.0,
        g    = 9.81,
        cp_d = (8.3144598 / 28.97e-3) / (2/7),
        cv_d = (8.3144598 / 28.97e-3) / (2/7) - 8.3144598 / 28.97e-3,
        xc   = 5000,
        yc   = 1000,
        zc   = 2000,
        rc   = 2000,
        xmax = 10000,
        ymax = 500,
        zmax = 10000,
        θₐ   = 2.0,
        cₛ   = 340,
        λ_coupler = ( 500 / 60 / 86400 ) #(L_airsea / τ_airsea)
    )


########
# Set up inital condition
########
# r(p, x, z)      = sqrt((x - p.xc)^2 + (z - p.zc)^2)
# Δθ(p, x, y, z)  = (r(p, x, z) < p.rc) ? ( p.θₐ * (1.0 - r(p, x, z) / p.rc) ) : 0.0
θ₀A(p, x, y, z)  = 250.0
π_exnerA(p, x, y, z)   = 1.0 - p.g / (p.cp_d * θ₀A(p, x, y, z) ) * z  

ρ₀A(p, x, y, z)  = p.pₒ / (p.R_d * θ₀A(p, x, y, z)) * (π_exnerA(p, x, y, z))^(p.cv_d/p.R_d)
# ρ₀(p, x, y, z)  = 1.0
ρθ₀A(p, x, y, z)     = ρ₀A(p, x, y, z) * θ₀A(p, x, y, z) 
ρu⃗₀A(p, x, y, z)     = @SVector [0,0,0]



# r(p, x, z)      = sqrt((x - p.xc)^2 + (z - p.zc)^2)
# Δθ(p, x, y, z)  = (r(p, x, z) < p.rc) ? ( p.θₐ * (1.0 - r(p, x, z) / p.rc) ) : 0.0
θ₀B(p, x, y, z)  = 300.0 + 10.0 * ( z /p.zmax ) 
π_exnerB(p, x, y, z)   = 1.0 - p.g / (p.cp_d * θ₀B(p, x, y, z) ) * z  

ρ₀B(p, x, y, z)  = p.pₒ / (p.R_d * θ₀B(p, x, y, z)) * (π_exnerB(p, x, y, z))^(p.cv_d/p.R_d)
ρθ₀B(p, x, y, z)     = ρ₀B(p, x, y, z) * θ₀B(p, x, y, z) 
ρu⃗₀B(p, x, y, z)     = @SVector [0,0,0]