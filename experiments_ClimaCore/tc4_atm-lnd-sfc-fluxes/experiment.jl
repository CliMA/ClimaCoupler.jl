# Minimal Example for SurfaceFluxes

# add https://github.com/CliMA/ClimaCore.jl/#main
# add https://github.com/CliMA/ClimaAtmos.jl/#main

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

import SurfaceFluxes
const SF = SurfaceFluxes

include("dummy_surface_fluxes.jl") 

parameters_ = (
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

########
# Experiment TC4
########
function exp_tc4(parameters::NamedTuple)
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

    ########
    # Set up rhs! and BCs
    ########

    # define model equations:
    include("atmos_rhs.jl")

    function ∑tendencies_lnd!(dT_sfc, T_sfc, (parameters, F_sfc), t)
        """
        Slab ocean:
        ∂_t T_sfc = F_sfc + G
        """
        p = parameters
        G = 0.0 # place holder for soil dynamics

        @. dT_sfc = (F_sfc[1] * p.C_p + G) / (p.h_s * p.ρ_s * p.c_s)
    end

    ########
    # Initiate prognostic variables
    ########

    # atmos IC state
    Yc = init_ekman_column_1d_c.(z_centers, Ref(parameters))
    Yf = init_ekman_column_1d_f.(z_faces, Ref(parameters))
    Y_atm_0 = ArrayPartition((Yc, Yf)) # (Y_centers, Y_faces, surface_fluxes)

    # land IC state
    Y_lnd_0 = [parameters.T_surf_lnd] # initiates lnd progostic var

    ics = (;
            atm = Y_atm_0,
            lnd = Y_lnd_0
            )

    ########
    # Set up time steppers
    ########

    # specify timestepping info
    stepping = (;
            Δt_min = parameters.Δt_min,
            timerange = parameters.timerange,
            Δt_cpl = max(parameters.nsteps_atm, parameters.nsteps_lnd) * parameters.Δt_min, # period of coupling cycle [s]
            odesolver = parameters.odesolver,
            nsteps_atm = parameters.nsteps_atm,
            nsteps_lnd = parameters.nsteps_lnd,
            )

    # Solve the ODE operator
    function coupler_solve!(stepping, ics, parameters)
        t = 0.0
        Δt_min  = stepping.Δt_min
        Δt_cpl  = stepping.Δt_cpl
        t_start = stepping.timerange[1]
        t_end   = stepping.timerange[2]
        nsteps_atm = stepping.nsteps_atm
        nsteps_lnd = stepping.nsteps_lnd

        # atmos copies of input vars (to be provided by coupler)
        atm_T_sfc_0 = ics.lnd
        atm_F_sfc_0 = [0.0, 0.0, 0.0]

        # SETUP ATMOS
        # put all prognostic variable arrays into a vector and ensure that solve can partition them
        Y_atm = ArrayPartition(( ics.atm.x[1], ics.atm.x[2] , atm_F_sfc_0))
        prob_atm = ODEProblem(∑tendencies_atm!, Y_atm, (t_start, t_end), (parameters, atm_T_sfc_0))
        integ_atm = init(
                            prob_atm,
                            stepping.odesolver,
                            dt = Δt_cpl / nsteps_atm,
                            saveat = parameters.saveat,)

        # land copies of coupler variables
        lnd_T_sfc = ics.lnd
        lnd_F_sfc = deepcopy(integ_atm.u.x[3] / Δt_cpl)
        
        # SETUP LAND
        prob_lnd = ODEProblem(∑tendencies_lnd!, lnd_T_sfc, (t_start, t_end), (parameters, lnd_F_sfc))
        integ_lnd = init(
                            prob_lnd,
                            stepping.odesolver,
                            dt = Δt_cpl / nsteps_lnd,
                            saveat = parameters.saveat,)

        # coupler stepping
        for t in (t_start : Δt_cpl : t_end)

            ## Atmos
            # pre_atmos
            integ_atm.u.x[3] .= [0.0, 0.0, 0.0] # surface flux to be accumulated

            # run atmos
            # NOTE: use (t - integ_atm.t) here instead of Δt_cpl to avoid accumulating roundoff error in our timestepping.
            step!(integ_atm, t - integ_atm.t, true)

            ## Land
            # pre_land
            integ_lnd.p[2] .= integ_atm.u.x[3] / Δt_cpl
            
            # run land
            step!(integ_lnd, t - integ_lnd.t, true)

            # post land
            integ_atm.p[2] .= deepcopy(integ_lnd.u) # get lnd_T_sfc from lnd and setup as aux variable for atm
        end

        return integ_atm, integ_lnd
    end

    ########
    # Run simulation
    ########

    integ_atm, integ_lnd = coupler_solve!(stepping, ics, parameters)
    sol_atm, sol_lnd = integ_atm.sol, integ_lnd.sol

    return sol_atm, sol_lnd
end

########
# Visualisation
########

function visualize(sol_atm, sol_lnd )
    ENV["GKSwstype"] = "nul"
    
    Plots.GRBackend()

    dirname = "heat"
    path = joinpath(@__DIR__, "output", dirname)
    mkpath(path)

    # animation
    # anim = Plots.@animate for u in sol_atm.u
    #     Plots.plot(u.x[1], xlim=(220,280))
    # end
    # Plots.mp4(anim, joinpath(path, "heat.mp4"), fps = 10)

    # atmos vertical profile at t=0 and t=end
    t0_ρθ = parent(sol_atm.u[1].x[1])[:,4]
    tend_ρθ = parent(sol_atm.u[end].x[1])[:,4]
    t0_u = parent(sol_atm.u[1].x[1])[:,2]
    tend_u = parent(sol_atm.u[end].x[1])[:,2]
    t0_v = parent(sol_atm.u[1].x[1])[:,3]
    tend_v = parent(sol_atm.u[end].x[1])[:,3]
    z_centers =  collect(1:1:length(tend_u))#parent(Fields.coordinate_field(center_space_atm))[:,1]
    Plots.png(Plots.plot([t0_ρθ tend_ρθ],z_centers, labels = ["t=0" "t=end"]), joinpath(path, "T_atm_height.png"))
    Plots.png(Plots.plot([t0_u tend_u],z_centers, labels = ["t=0" "t=end"]), joinpath(path, "u_atm_height.png"))
    Plots.png(Plots.plot([t0_v tend_v],z_centers, labels = ["t=0" "t=end"]), joinpath(path, "v_atm_height.png"))

    # time evolution (convert to total_energy in final test when operators fixed)
    atm_sum_u_t = [sum(parent(u.x[1])[:,4]) for u in sol_atm.u] ./ parameters_.atmos_Nz .* parameters_.atmos_Lz * parameters_.C_p # J / m2
    lnd_sfc_u_t = [u[1] for u in sol_lnd.u] * parameters_.h_s * parameters_.ρ_s * parameters_.c_s # J / m2

    v1 = lnd_sfc_u_t .- lnd_sfc_u_t[1]
    v2 = atm_sum_u_t .- atm_sum_u_t[1]
    Plots.png(Plots.plot(sol_lnd.t, [v1 v2 v1+v2], labels = ["lnd" "atm" "tot"]), joinpath(path, "energy_both_surface_time.png"))

    # relative error
    total = atm_sum_u_t + lnd_sfc_u_t
    rel_error = (total .- total[1]) / mean(total)
    Plots.png(Plots.plot(sol_lnd.t, rel_error, labels = ["tot"]), joinpath(path, "rel_error_surface_time.png"))

    function linkfig(figpath, alt = "")
        # buildkite-agent upload figpath
        # link figure in logs if we are running on CI
        if get(ENV, "BUILDKITE", "") == "true"
            artifact_url = "artifact://$figpath"
            print("\033]1338;url='$(artifact_url)';alt='$(alt)'\a\n")
        end
    end

    linkfig("output/$(dirname)/TC4_end.png", "TC4 End Simulation")
end


sol_atm, sol_lnd = exp_tc4(parameters_)
#= 
using Statistics
import Plots
visualize(sol_atm, sol_lnd );
=#
