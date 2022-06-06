const FT = Float64
import SciMLBase: step!
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!
import ClimaCore.Utilities: PlusHalf
using DiffEqCallbacks

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

@info "Init Models and Maps"

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

atm_boundary = Spaces.level(atm_domain.hv_face_space, PlusHalf(0))

# Build remapping operators
atm_to_ocn = Operators.LinearRemap(ocn_domain, atm_boundary)
atm_to_lnd = Operators.LinearRemap(lnd_domain, atm_boundary)
ocn_to_atm = Operators.LinearRemap(atm_boundary, ocn_domain)
lnd_to_atm = Operators.LinearRemap(atm_boundary, lnd_domain)

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

function cpl_run()

    # timestepping
    t_start, t_end = (0.0, 5e2)
    cpl_Δt = 0.1
    saveat = 1e2
    atm_nsteps, ocn_nsteps, lnd_nsteps = (5, 1, 1)

    # initialize
    atm_T_sfc = Operators.remap(ocn_to_atm, ocn_Y_default.T_sfc) .+ Operators.remap(lnd_to_atm, lnd_Y_default.T_sfc) # masked arrays; regrid to atm grid
    atm_F_sfc = Fields.zeros(atm_boundary)

    atm_Y = Fields.FieldVector(Yc = atm_Y_default.Yc, ρw = atm_Y_default.ρw, F_sfc = atm_F_sfc)
    atm_prob = ODEProblem(
        atm_rhs!,
        atm_Y,
        (t_start, t_end),
        (cpl_p = cpl_parameters, T_sfc = atm_T_sfc, bc = atm_bc, domain = atm_domain),
    )
    atm_integ = init(
        atm_prob,
        SSPRK33(),
        dt = cpl_Δt / atm_nsteps,
        saveat = saveat,
        progress = true,
        progress_message = (dt, u, params, t) -> t,
        callback = dss_callback,
    )

    ocn_F_sfc = Fields.zeros(ocn_domain)

    ocn_Y = Fields.FieldVector(T_sfc = ocn_Y_default.T_sfc)
    ocn_prob = ODEProblem(ocn_rhs!, ocn_Y, (t_start, t_end), (cpl_parameters, F_sfc = ocn_F_sfc))
    ocn_integ = init(ocn_prob, SSPRK33(), dt = cpl_Δt / ocn_nsteps, saveat = saveat)

    lnd_F_sfc = Fields.zeros(lnd_domain)

    lnd_Y = Fields.FieldVector(T_sfc = lnd_Y_default.T_sfc)
    lnd_prob = ODEProblem(lnd_rhs!, lnd_Y, (t_start, t_end), (cpl_parameters, F_sfc = lnd_F_sfc))
    lnd_integ = init(lnd_prob, SSPRK33(), dt = cpl_Δt / lnd_nsteps, saveat = saveat)

    @info "Run model"
    # coupler stepping
    for t in (t_start:cpl_Δt:t_end)

        ## Atmos
        # pre: reset flux accumulator
        atm_F_sfc = atm_integ.u.F_sfc
        atm_F_sfc .= atm_F_sfc .* 0.0    # reset surface flux to be accumulated

        # run 
        # NOTE: use (t - integ_atm.t) here instead of Δt_cpl to avoid accumulating roundoff error in our timestepping.
        step!(atm_integ, t - atm_integ.t, true)

        ## Ocean
        # pre: get accumulated flux from atmos
        ocn_F_sfc = ocn_integ.p.F_sfc
        Operators.remap!(ocn_F_sfc, atm_to_ocn, atm_integ.u.F_sfc ./ cpl_Δt)

        # run
        step!(ocn_integ, t - ocn_integ.t, true)
        # post: send ocean surface temp to atmos
        Operators.remap!(atm_integ.p.T_sfc, ocn_to_atm, ocn_integ.u.T_sfc)

        ## Land
        # pre: get accumulated flux from atmos
        lnd_F_sfc = lnd_integ.p.F_sfc
        Operators.remap!(lnd_F_sfc, atm_to_lnd, atm_integ.u.F_sfc ./ cpl_Δt)

        # run
        step!(lnd_integ, t - lnd_integ.t, true)
        # post: send ocean surface temp to atmos
        atm_integ.p.T_sfc .+= Operators.remap(lnd_to_atm, lnd_integ.u.T_sfc)
    end
    @info "Simulation Complete"
    return atm_integ, ocn_integ, lnd_integ
end

atm_integ, ocn_integ, lnd_integ = cpl_run()

sol = atm_integ.sol

# dirname = "sea_breeze_2d"
# path = joinpath(@__DIR__, "output", dirname)
# mkpath(path)

# using JLD2
# save(joinpath(path, "last_sim.jld2"), "atm_integ", atm_integ, "ocn_integ", ocn_integ, "lnd_integ", lnd_integ)

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
