# Minimal Example for SurfaceFluxes
using ClimaCore.Geometry

using ClimaCore: DataLayouts, Operators, Geometry
using ClimaAtmos.Simulations: Simulation
using ClimaCore: Fields, Domains, Topologies, Meshes, Spaces

using IntervalSets
using UnPack
using OrdinaryDiffEq: SSPRK33
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33
using Logging: global_logger
using TerminalLoggers: TerminalLogger

using RecursiveArrayTools
using LinearAlgebra
using OrdinaryDiffEq
using Random

using DocStringExtensions

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const CLIMAparam_set = EarthParameterSet()
import CLIMAParameters

abstract type AbstractNonlinearSolverMethod{FT} end
abstract type  AbstractTolerance{FT} end

const FTypes = Union{AbstractFloat, AbstractArray}

global_logger(TerminalLogger())

const CI = !isnothing(get(ENV, "CI", nothing))

import SciMLBase: step!
using Printf

const FT = Float64

# define model equations:
include("atmos_rhs_minimal_example.jl")
include("SurfaceFluxes.jl")
include("dummy_surface_fluxes.jl") 

parameters = (
        # atmos domain
        atmos_Lz = 10,
        atmos_Nz = 50,
        # atmos physics (overwriting CliMAParameters defaults?)
        MSLP = 1e5, # mean sea level pressure [Pa]
        grav = FT(9.8), # gravitational acceleration [m /s^2]
        R_d = 287.058, # R dry (gas constant / mol mass dry air)  [J / K / kg]
        C_p = 287.058 * 1.4 / (1.4 - 1), # heat capacity at constant pressure [J / K / kg]
        C_v = 287.058 / (1.4 - 1), # heat capacity at constant volume [J / K / kg]
        R_m = 87.058, # moist R, assumed to be dry [J / K / kg]
        f = 7.29e-5, # Coriolis parameters [1/s]
        ν = 0.1, #0.01 # viscosity, diffusivity
        Ch = 0.0015, # bulk transfer coefficient for sensible heat
        Cd = 0.01 / (2e2 / 30.0), #drag coeff
        ug = 1.0,
        vg = 0.0,
        d = sqrt(2.0 * 0.01 / 5e-5), #?
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
        T_surf = 260.0,
        T_min_ref = 230.0,
        u0 = 1.0,
        v0 = 0.0,
        w0 = 0.0,
    )

# Get surface temperature from the land simulation

start_time = 0.0
stop_time = 1.0
Δt_min  = 0.02

########
# Set up atmos domain
########

domain_atm  = Domains.IntervalDomain(0, parameters.atmos_Lz, x3boundary = (:bottom, :top)) # struct
mesh_atm = Meshes.IntervalMesh(domain_atm, nelems = parameters.atmos_Nz) # struct, allocates face boundaries to 5,6
center_space_atm = Spaces.CenterFiniteDifferenceSpace(mesh_atm) # collection of the above, discretises space into FD and provides coords
face_space_atm = Spaces.FaceFiniteDifferenceSpace(center_space_atm)

# Initialize the atmospheric states Yc and Yf
z_centers = Fields.coordinate_field(center_space_atm)
z_faces = Fields.coordinate_field(face_space_atm)

Yc = init_centers.(z_centers, Ref(parameters))
Yf = init_faces.(z_faces, Ref(parameters))

# Put all prognostic variable arrays into a vector and ensure that solve can partition them
Y_atm = ArrayPartition((Yc, Yf, zeros(3)))


prob_atm = ODEProblem(∑tendencies_atm!, Y_atm, (start_time, stop_time), (parameters, [parameters.T_surf]))
simulation = init(prob_atm, SSPRK33(), dt = Δt_min, saveat = 1 * Δt_min)

function coupler_solve()
    for t in (start_time : Δt_min : stop_time)
        step!(simulation, t - simulation.t, true)
    end
    return simulation
end
# while coupled_sim.clock.time < stop_time
#     step!(simulation, stop_time - atmos_sim.t, true)
# end
      



