const FT = Float64
import SciMLBase: step!
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!
import ClimaCore.Utilities: PlusHalf
using DiffEqCallbacks

push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", ".."))
using ClimaCoupler

include("atmos_rhs.jl")
include("ocean_rhs.jl")
include("land_rhs.jl")
include("coupledbc.jl")

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
    # atmos parameters
    atm_μ = FT(0.0001), # diffusion coefficient
    atm_T_top = FT(280.0), # fixed temperature at the top of the domain_atm
    atm_T_ini = atm_T_ini, # initial condition of at temperature (isothermal) [K]
    MSLP = MSLP, # mean sea level pressure
    grav = grav, # gravitational constant
    R_d = R_d, # R dry (gas constant / mol mass dry air)
    γ = γ, # heat capacity ratio
    C_p = C_p, # heat capacity at constant pressure
    C_v = C_v, # heat capacity at constant volume
    R_m = R_m, # moist R, assumed to be dry
    # land slab parameters
    lnd_h = FT(0.5), # depth of slab layer [m]
    lnd_ρ = FT(1500), # density [kg m^-3]
    lnd_c = FT(800), # specific heat [J K^-1 kg^-1]
    lnd_T_ini = FT(260.0), # initial condition of at temperature (isothermal) [K]
    # ocean slab parameters
    ocn_h = FT(0.5), # depth of slab layer [m]
    ocn_ρ = FT(1025), # density [kg m^-3]
    ocn_c = FT(3850), # specific heat [J K^-1 kg^-1]
    ocn_T_ini = FT(260.0), # initial condition of at temperature (isothermal) [K]
    # coupling parameters
    C_H = FT(0.0015),
)

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

# timestepping
t_start, t_end = (0.0, 1.0)
Δt_coupled = 0.1
saveat = 1e2
atm_nsteps, ocn_nsteps, lnd_nsteps = (5, 1, 1)

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

# Build remapping operators
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
    Operators.remap(maps.land_to_atmos, lnd_Y_default.T_sfc) # masked arrays; regrid to atm grid
atm_F_sfc = Fields.zeros(atm_boundary)
ocn_F_sfc = Fields.zeros(ocn_domain)
lnd_F_sfc = Fields.zeros(lnd_domain)

# init models
atm_Y = Fields.FieldVector(Yc = atm_Y_default.Yc, ρw = atm_Y_default.ρw, F_sfc = atm_F_sfc)
atm_p = (cpl_p = cpl_parameters, T_sfc = atm_T_sfc, bc = atm_bc)
atmos = AtmosSimulation(atm_Y, t_start, Δt_coupled / atm_nsteps, t_end, SSPRK33(), atm_p, saveat, dss_callback)

ocn_Y = Fields.FieldVector(T_sfc = ocn_Y_default.T_sfc)
ocn_p = (cpl_parameters, F_sfc = ocn_F_sfc)
ocean = OceanSimulation(ocn_Y, t_start, Δt_coupled / ocn_nsteps, t_end, SSPRK33(), ocn_p, saveat)

lnd_Y = Fields.FieldVector(T_sfc = lnd_Y_default.T_sfc)
lnd_p = (cpl_parameters, F_sfc = lnd_F_sfc)
land = LandSimulation(lnd_Y, t_start, Δt_coupled / lnd_nsteps, t_end, SSPRK33(), lnd_p, saveat)

# coupled simulation
struct AOLCoupledSimulation{
    A <: AtmosSimulation,
    O <: OceanSimulation,
    L <: LandSimulation,
    C <: ClimaCoupler.CouplerState,
} <: ClimaCoupler.AbstractCoupledSimulation
    # Atmosphere Simulation
    atmos::A
    # Ocean Simulation
    ocean::O
    # Land Simulation
    land::L
    # Coupler storage
    coupler::C
end

# init coupler fields and maps
coupler = CouplerState(Δt_coupled)
coupler_add_field!(coupler, :T_sfc_ocean, ocean.integrator.u.T_sfc; write_sim = ocean)
coupler_add_field!(coupler, :T_sfc_land, land.integrator.u.T_sfc; write_sim = land)
coupler_add_field!(coupler, :F_sfc, atmos.integrator.u.F_sfc; write_sim = atmos)
for (name, map) in pairs(maps)
    coupler_add_map!(coupler, name, map)
end

sim = AOLCoupledSimulation(atmos, ocean, land, coupler)

# step for sims built on OrdinaryDiffEq
# NOTE: use (t - integ_atm.t) here instead of Δt_cpl to avoid accumulating roundoff error in our timestepping.
function step!(sim::ClimaCoupler.AbstractSimulation, t_stop)
    Δt = t_stop - sim.integrator.t
    step!(sim.integrator, Δt, true)
end

function cpl_run(simulation::AOLCoupledSimulation)
    @info "Run model"
    @unpack atmos, ocean, land, coupler = simulation
    Δt_coupled = coupler.Δt_coupled
    # coupler stepping
    for t in ((t_start + Δt_coupled):Δt_coupled:t_end)
        ## Atmos
        coupler_pull!(atmos, coupler)
        step!(atmos, t)
        coupler_push!(coupler, atmos)

        ## Ocean
        coupler_pull!(ocean, coupler)
        step!(ocean, t)
        coupler_push!(coupler, ocean)

        ## Land
        coupler_pull!(land, coupler)
        step!(land, t)
        coupler_push!(coupler, land)
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
