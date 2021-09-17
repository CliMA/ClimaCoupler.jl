# Experiments

using Test

parameters = (
        # timestepping parameters 
        Δt_min = 0.01, # minimum model timestep [s]
        timerange = (0.0, 10.0),  # start time and end time [s]
        odesolver = SSPRK33(), # timestepping method from DifferentialEquations.jl (used in both models here)
        nsteps_atm = 1, # no. time steps of atm before coupling 
        nsteps_lnd = 1, # no. time steps of lnd before coupling 
        saveat = 0.2, # interval at which to save diagnostics [s]

        # atmos domain
        atmos_Lz = 300,
        atmos_Nz = 10,
        
        # atmos physics (overwriting CliMAParameters-v2 defaults?)
        MSLP = FT(1e5), # mean sea level pressure
        grav = FT(9.8), # gravitational constant
        R_d = FT(287.058), # R dry (gas constant / mol mass dry air)
        C_p = FT(287.058 * 7 / 2), # heat capacity at constant pressure
        C_v = FT(287.058 * 5 / 2), # heat capacity at constant volume
        R_m = FT(287.058), # moist R, assumed to be dry
        f = FT(5e-5), # Coriolis parameters
        ν = FT(0.01),
        uvg = Geometry.Cartesian12Vector(FT(1.0), FT(0.0)),
        T_min_ref = FT(230.0),
        T_surf_atm = FT(270.0),
        u0 = FT(1.0),
        v0 = FT(0.0),
        w0 = FT(0.0),

        # soil slab
        h_s  = 100,     # depth of the modelled soil model [m]
        c_s  = 800,     # specific heat for land (soil)  [J / K / kg]
        κ_s  = 0.0,     # soil conductivity [W / m / K] (set to 0 for energy conservation checks)
        T_h  = 280,     # temperature of soil at depth h [K]
        ρ_s  = 1500,    # density for land (soil) [kg / m^3]
        T_surf_lnd = 280.0, # initial surface temperature [K]

        # radiation parameters for DryBulkFormulaWithRadiation() SurfaceFluxType
        τ    = 0.9,     # atmospheric transmissivity
        α    = 0.5,     # surface albedo
        σ    = 5.67e-8, # Steffan Boltzmann constant [kg / s^3 / K^4]
        g_a  = 0.06,    # aerodynamic conductance for heat transfer [kg / m^2 / s]
        ϵ    = 0.98,    # broadband emissivity / absorptivity of the surface
        F_a  = 0.0,     # downward LW flux from the atmosphere [W / m^2]
        F_sol = 1361,   # incoming solar TOA radiation [W / m^2]
        τ_d   = 10,     # idealized daily cycle period [s]

        # surface fluxes
        λ = FT(0.01),#FT(1e-5)    # coupling transfer coefficient for LinearRelaxation() SurfaceFluxType 
        Ch = 0.0015, # bulk transfer coefficient for sensible heat
        Cd = 0.0015, # drag coefficient
    )

@testset "Coupler Interface" begin

    PWD = pwd()
    Pkg.activate(PWD*"/experiments_ClimaCore/tc4_atm-lnd-sfc-fluxes/")
    include(PWD*"/experiments_ClimaCore/tc4_atm-lnd-sfc-fluxes/experiment.jl")

    # check if runs
    sol_atm, sol_lnd = exp_tc4(parameters)

    # conservation checks - add when div operator fixed

end
