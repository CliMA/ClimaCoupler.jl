# # Coupled Sea Breeze Driver

# In this tutorial we demonstrate the coupling of three component models
# (atmosphere, ocean, and land) to drive the sea breeze. The primary parts
# of the ClimaCoupler interface are used and discussed.

const FT = Float64
import SciMLBase: step!
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!
import ClimaCore.Utilities: PlusHalf
using DiffEqCallbacks

# enable broadcasting with mismatched spaces
import ClimaCore: Fields, Operators
Fields.allow_mismatched_diagonalized_spaces() = true
Operators.allow_mismatched_fd_spaces() = true

push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", ".."))
using ClimaCoupler

#=
## Component Models
Coupled simulations can re-use tendency methods developed for standalone simulations.
This is achieved by multiple dispatch, where methods that deal with boundaries
dispatch off of a coupled boundary type. This minimizes the necessary code that
must be specialized for a coupled run as only special boundary conditions must
be written. Here, the atmosphere has special boundary conditions for coupling,
while the ocean and land tendencies are unaltered.
=#

include("atmos_rhs.jl")
include("ocean_rhs.jl")
include("land_rhs.jl")

# model parameters that the coupler overwrites 
const atm_T_ini = FT(270.0)
const MSLP = FT(1e5)
const grav = FT(9.8)
const R_d = FT(287.058)
const γ = FT(1.4)
const C_p = FT(R_d * γ / (γ - 1))
const C_v = FT(R_d / (γ - 1))
const R_m = R_d
cpl_parameters = (
    ## atmos parameters
    atm_μ = FT(0.0001), ## diffusion coefficient
    atm_T_top = FT(280.0), ## fixed temperature at the top of the domain_atm
    atm_T_ini = atm_T_ini, ## initial condition of at temperature (isothermal) [K]
    MSLP = MSLP, ## mean sea level pressure
    grav = grav, ## gravitational constant
    R_d = R_d, ## R dry (gas constant / mol mass dry air)
    γ = γ, ## heat capacity ratio
    C_p = C_p, ## heat capacity at constant pressure
    C_v = C_v, ## heat capacity at constant volume
    R_m = R_m, ## moist R, assumed to be dry
    ## land slab parameters
    lnd_h = FT(0.5), ## depth of slab layer [m]
    lnd_ρ = FT(1500), ## density [kg m^-3]
    lnd_c = FT(800), ## specific heat [J K^-1 kg^-1]
    lnd_T_ini = FT(260.0), ## initial condition of at temperature (isothermal) [K]
    ## ocean slab parameters
    ocn_h = FT(0.5), ## depth of slab layer [m]
    ocn_ρ = FT(1025), ## density [kg m^-3]
    ocn_c = FT(3850), ## specific heat [J K^-1 kg^-1]
    ocn_T_ini = FT(260.0), ## initial condition of at temperature (isothermal) [K]
    ## coupling parameters
    C_H = FT(0.0015),
)

# DSS callback
function make_dss_func()
    function _dss!(x::Fields.Field)
        Spaces.weighted_dss!(x)
    end
    function _dss!(::Any)
        nothing
    end
    dss_func(Y, t, integrator) = foreach(_dss!, Fields._values(Y))
    return dss_func
end
dss_func = make_dss_func()
dss_callback = FunctionCallingCallback(dss_func, func_start = true)

@info "Init Models and Maps"

# ## Timestepping
# The coupled simulation synchronizes the component models at a coupling time step,
# `Δt_cpl`. Within that step, components may substep - each component specifies a
# number of substeps to take within `Δt_cpl`: `atm_nsteps, ocn_nsteps, lnd_nsteps`.
t_start, t_end = (0.0, 1.0)
Δt_cpl = 0.1
saveat = 1e2
atm_nsteps, ocn_nsteps, lnd_nsteps = (5, 1, 1)

# Initialize Fields
atm_Y_default, atm_bc, atm_domain = atm_init(
    xmin = -500,
    xmax = 500,
    zmin = 0,
    zmax = 1000,
    npoly = 4,
    helem = 20,
    velem = 20,
    bc = (ρθ = (bottom = CoupledFlux(), top = ZeroFlux()),),
)

ocn_Y_default, ocn_domain = ocn_init(xmin = -500, xmax = 0, helem = 10, npoly = 0)

lnd_Y_default, lnd_domain = lnd_init(xmin = 0, xmax = 500, helem = 10, npoly = 0)

# ## Remapping
# Because models may live on different grids, remapping is necessary at the boundaries.
# Maps between coupled components must be constructed for each interacting pair.
atm_boundary = Spaces.level(atm_domain.hv_face_space, PlusHalf(0))

maps = (
    atmos_to_ocean = Operators.LinearRemap(ocn_domain, atm_boundary),
    atmos_to_land = Operators.LinearRemap(lnd_domain, atm_boundary),
    ocean_to_atmos = Operators.LinearRemap(atm_boundary, ocn_domain),
    land_to_atmos = Operators.LinearRemap(atm_boundary, lnd_domain),
)

# initialize coupling fields
atm_T_sfc =
    Operators.remap(maps.ocean_to_atmos, ocn_Y_default.T_sfc) .+
    Operators.remap(maps.land_to_atmos, lnd_Y_default.T_sfc) ## masked arrays; regrid to atm grid
atm_F_sfc = Fields.zeros(atm_boundary)
ocn_F_sfc = Fields.zeros(ocn_domain)
lnd_F_sfc = Fields.zeros(lnd_domain)

# ## Simulations
# Each component is wrapped as a Simulation, which contains both the model (tendency)
# and the time-stepping information (solver, step size, etc). Simulations are the standard
# structures that the coupler works with, enabling dispatch. Here, we create three simulations:
# `AtmosSimulation`, `OceanSimulation`, and `LandSimulation`.
atm_Y = Fields.FieldVector(Yc = atm_Y_default.Yc, ρw = atm_Y_default.ρw, F_sfc = atm_F_sfc)
atm_p = (cpl_p = cpl_parameters, T_sfc = atm_T_sfc, bc = atm_bc)
atmos = AtmosSimulation(atm_Y, t_start, Δt_cpl / atm_nsteps, t_end, SSPRK33(), atm_p, saveat, dss_callback)

ocn_Y = Fields.FieldVector(T_sfc = ocn_Y_default.T_sfc)
ocn_p = (cpl_parameters, F_sfc = ocn_F_sfc)
ocean = OceanSimulation(ocn_Y, t_start, Δt_cpl / ocn_nsteps, t_end, SSPRK33(), ocn_p, saveat)

lnd_Y = Fields.FieldVector(T_sfc = lnd_Y_default.T_sfc)
lnd_p = (cpl_parameters, F_sfc = lnd_F_sfc)
land = LandSimulation(lnd_Y, t_start, Δt_cpl / lnd_nsteps, t_end, SSPRK33(), lnd_p, saveat)

# Additionally, we create a coupled simulation that contains the component simulations
# and the coupled time-stepping information.
struct AOLCoupledSimulation{
    FT,
    A <: AtmosSimulation,
    O <: OceanSimulation,
    L <: LandSimulation,
    C <: ClimaCoupler.CouplerState,
} <: ClimaCoupler.AbstractCoupledSimulation
    ## Atmosphere Simulation
    atmos::A
    ## Ocean Simulation
    ocean::O
    ## Land Simulation
    land::L
    ## Coupler storage
    coupler::C
    ## The coupled time step size
    Δt::FT
end

# ## The Coupler
# The `CouplerState` is a coupling struct used to store pointers or copies of the
# shared boundary information. All components are coupled by updating or accessing
# data in this `CouplerState`; component models do not directly interface with one another,
# only through the coupler.
coupler = CouplerState()
coupler_add_field!(coupler, :T_sfc_ocean, ocean.integrator.u.T_sfc; write_sim = ocean)
coupler_add_field!(coupler, :T_sfc_land, land.integrator.u.T_sfc; write_sim = land)
coupler_add_field!(coupler, :F_sfc, atmos.integrator.u.F_sfc; write_sim = atmos)
for (name, map) in pairs(maps)
    coupler_add_map!(coupler, name, map)
end

sim = AOLCoupledSimulation(atmos, ocean, land, coupler, Δt_cpl)

# step for sims built on OrdinaryDiffEq
function step!(sim::ClimaCoupler.AbstractSimulation, t_stop)
    Δt = t_stop - sim.integrator.t
    step!(sim.integrator, Δt, true)
end

function cpl_run(simulation::AOLCoupledSimulation)
    @info "Run model"
    @unpack atmos, ocean, land, Δt = simulation
    Δt_cpl = Δt
    ## coupler stepping
    for t in ((t_start + Δt_cpl):Δt_cpl:t_end)

        ## Atmos
        ## pre: reset flux accumulator
        atmos.integrator.u.F_sfc .= 0.0 # reset surface flux to be accumulated
        # don't want to alloc here..
        T_sfc_ocean = coupler_get(coupler, :T_sfc_ocean, atmos)
        T_sfc_land = coupler_get(coupler, :T_sfc_land, atmos)
        atmos.integrator.p.T_sfc .= T_sfc_land .+ T_sfc_ocean

        # run
        # NOTE: use (t - integ_atm.t) here instead of Δt_cpl to avoid accumulating roundoff error in our timestepping.
        step!(atmos, t)
        coupler_put!(coupler, :F_sfc, atmos.integrator.u.F_sfc, atmos)

        ## Ocean
        # pre: get accumulated flux from atmos
        ocn_F_sfc = ocean.integrator.p.F_sfc
        ocn_F_sfc .= coupler_get(coupler, :F_sfc, ocean) ./ Δt_cpl

        # run
        step!(ocean, t)
        # post: send ocean surface temp to atmos
        coupler_put!(coupler, :T_sfc_ocean, ocean.integrator.u.T_sfc, ocean)

        ## Land
        # pre: get accumulated flux from atmos
        lnd_F_sfc = land.integrator.p.F_sfc
        lnd_F_sfc .= coupler_get(coupler, :F_sfc, land) ./ Δt_cpl

        # run
        step!(land, t)
        # post: send land surface temp to atmos
        coupler_put!(coupler, :T_sfc_land, land.integrator.u.T_sfc, land)
    end
    @info "Simulation Complete"
end

cpl_run(sim)

# sol = sim.atmos.integrator.sol

# dirname = "sea_breeze_2d"
# path = joinpath(@__DIR__, "output", dirname)
# mkpath(path)

# using JLD2
# save(joinpath(path, "last_sim.jld2"), "coupled_sim", sim)

# post-processing
# import Plots, ClimaCorePlots
# Plots.GRBackend()
# interp = 5

# anim = Plots.@animate for u in sol.u
#     Plots.contourf(u.Yc.ρθ ./ u.Yc.ρ)
# end
# Plots.mp4(anim, joinpath(path, "theta.mp4"), fps = 20)

#  If2c = Operators.InterpolateF2C()
#  anim = Plots.@animate for u in sol.u
#      Plots.contourf(If2c.(u.ρw) ./ u.Yc.ρ)
#  end
#  Plots.mp4(anim, joinpath(path, "vel_w.mp4"), fps = 20)

#  anim = Plots.@animate for u in sol.u
#      Plots.contourf(u.Yc.ρuₕ ./ u.Yc.ρ)
#  end
#  Plots.mp4(anim, joinpath(path, "vel_u.mp4"), fps = 20)
