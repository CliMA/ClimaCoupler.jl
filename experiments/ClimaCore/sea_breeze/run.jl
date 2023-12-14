# # Coupled Sea Breeze

#=
## Overview

This sea breeze simulation consists of an atmosphere above ocean and land thermal slabs.
The difference in heating between the land and ocean components drives circulation:
cool ocean air flows towards the land at the surface while warm air over land rises
and flows over the ocean.

In this tutorial we demonstrate the coupling of three component models
(atmosphere, ocean, and land) to drive the sea breeze. The primary parts
of the ClimaCoupler interface are used and discussed.
=#

# Load utilities for running coupled simulation
include("../CoupledSims/coupled_sim.jl")

using SciMLBase: ODEProblem, savevalues!, solve, init, CallbackSet #hide
import SciMLBase #hide
import ClimaTimeSteppers as CTS #hide
import ClimaCore.Utilities: PlusHalf #hide
import ClimaCore.Spaces as Spaces
using DiffEqCallbacks #hide

## enable broadcasting with mismatched spaces #hide
import ClimaCore: Fields, Operators #hide
Fields.allow_mismatched_diagonalized_spaces() = true #hide
Operators.allow_mismatched_fd_spaces() = true #hide
#hide
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "..")) #hide
using ClimaCoupler #hide
#hide
const FT = Float64 #hide

#=
## Model Initialization
### Component Models
Component models are the building blocks of coupled models. They are often developed
independently from one another and can be executed by themselves as "standalone" simulations.
The coupler is used to combine these components into coupled simulations.
Importantly, coupled simulations can re-use tendency methods developed for standalone simulations,
maximizing code reuse and minimizing the necessary code that
must be specialized for a coupled run--only special boundary conditions must
be written. This is achieved by multiple dispatch, where methods that deal with boundaries
dispatch off of a coupled boundary type. Here, the atmosphere has special boundary conditions
for coupling while the ocean and land tendencies are unaltered. See the atmospheric model page
for more details.

In a more mature CliMA ecosystem, the following include statements would be replaced
by `using` statements for the relevant component packages.
=#
include("atmos_rhs.jl")
include("ocean_rhs.jl")
include("land_rhs.jl")

## model parameters
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
    ## land slab parameters
    lnd_h = FT(0.5), # depth of slab layer [m]
    lnd_ρ = FT(1500), # density [kg m^-3]
    lnd_c = FT(800), # specific heat [J K^-1 kg^-1]
    lnd_T_ini = FT(260.0), # initial condition of at temperature (isothermal) [K]
    ## ocean slab parameters
    ocn_h = FT(0.5), # depth of slab layer [m]
    ocn_ρ = FT(1025), # density [kg m^-3]
    ocn_c = FT(3850), # specific heat [J K^-1 kg^-1]
    ocn_T_ini = FT(260.0), # initial condition of at temperature (isothermal) [K]
    ## coupling parameters
    C_H = FT(0.0015),
)

## DSS callback
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

#=
## Initialization
The coupled simulation synchronizes the component models at a coupling time step,
`Δt_cpl`. Within that step, components may substep - each component specifies a
number of substeps to take within `Δt_cpl`: `atm_nsteps, ocn_nsteps, lnd_nsteps`.

Component model states are initialized via the initialization methods each component
would use in standalone mode. These states will be modified to reflect the full coupled
system before executing the simulation.
=#
@info "Init Models and Maps"

t_start, t_end = (0.0, 1e4)
Δt_coupled = 0.1
saveat = 1.0
atm_nsteps, ocn_nsteps, lnd_nsteps = (5, 1, 1)

## Initialize Models
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

#=
## Remapping
Because models may live on different grids, remapping is necessary at the boundaries.
Maps between coupled components must be constructed for each interacting pair. Remapping
utilities are imported from `ClimaCore.Operators`.
=#
atm_boundary = Spaces.level(atm_domain.hv_face_space, PlusHalf(0))

maps = (
    atmos_to_ocean = Operators.LinearRemap(ocn_domain, atm_boundary),
    atmos_to_land = Operators.LinearRemap(lnd_domain, atm_boundary),
    ocean_to_atmos = Operators.LinearRemap(atm_boundary, ocn_domain),
    land_to_atmos = Operators.LinearRemap(atm_boundary, lnd_domain),
)

## initialize coupling fields
atm_T_sfc =
    Operators.remap(maps.ocean_to_atmos, ocn_Y_default.T_sfc) .+
    Operators.remap(maps.land_to_atmos, lnd_Y_default.T_sfc) # masked arrays; regrid to atm grid
atm_F_sfc = Fields.zeros(atm_boundary)
ocn_F_sfc = Fields.zeros(ocn_domain)
lnd_F_sfc = Fields.zeros(lnd_domain)

#=
## Simulations
Each component is wrapped as a `Sim`, which contains both the model (tendency)
and the time-stepping information (solver, step size, etc). Sims are the standard
structures that the coupler works with, enabling dispatch of coupler methods.
Here, we create three simulations: `AtmosSim`, `OceanSim`, and `LandSim`.
=#
atm_Y = Fields.FieldVector(Yc = atm_Y_default.Yc, ρw = atm_Y_default.ρw, F_sfc = atm_F_sfc)
atm_p = (cpl_p = cpl_parameters, T_sfc = atm_T_sfc, bc = atm_bc)
atmos = AtmosSim(atm_Y, t_start, Δt_coupled / atm_nsteps, t_end, CTS.RK4(), atm_p, saveat, dss_callback)

ocn_Y = Fields.FieldVector(T_sfc = ocn_Y_default.T_sfc)
ocn_p = (cpl_parameters, F_sfc = ocn_F_sfc)
ocean = OceanSim(ocn_Y, t_start, Δt_coupled / ocn_nsteps, t_end, CTS.RK4(), ocn_p, saveat)

lnd_Y = Fields.FieldVector(T_sfc = lnd_Y_default.T_sfc)
lnd_p = (cpl_parameters, F_sfc = lnd_F_sfc)
land = LandSim(lnd_Y, t_start, Δt_coupled / lnd_nsteps, t_end, CTS.RK4(), lnd_p, saveat)

# Additionally, we create a coupled simulation that contains the component simulations
# and the coupled time-stepping information.
struct AOLCoupledSim{A <: AtmosSim, O <: OceanSim, L <: LandSim, C <: CouplerState} <: AbstractCoupledSim
    ## Atmosphere Simulation
    atmos::A
    ## Ocean Simulation
    ocean::O
    ## Land Simulation
    land::L
    ## Coupler storage
    coupler::C
end

#=
`step!` is a key method within the Sims interface. It advances a simulation
to the specified `t_stop`, with that simulation advancing by its own internal step
size to reach the specified time. Each simulation type should specify its own step
method, allowing components to have different time integration backends. Here, all
components are using SciMLBase integrators and can share the same `step!` method.
=#
function step!(sim::AbstractSim, t_stop)
    Δt = t_stop - sim.integrator.t
    SciMLBase.step!(sim.integrator)
end

#=
## The Coupler
The `CouplerState` is a coupling struct used to store pointers or copies of the
shared boundary information. All components are coupled by updating or accessing
data in this `CouplerState`; component models do not directly interface with one another,
only through the coupler.

After creating the `CouplerState` object, coupled fields can be registered index the
coupler via the `coupler_add_field!` method. This field is then accessible by `coupler_get`
methods and can be updated via the `coupler_put!` methods.

Similarly, the `coupler_add_map!` method registers remapping operators in the coupler. To
provide automatic remapping, there is a strict name convention for remap operators: a map
from SimA to SimB (where `ClimaCoupler.name` returns `:simA` and `:simB`,
respectively) must be named `simA_to_simB` so that the correct operator can be used.

Here, the models are coupled through heat transfer at the surface. This heat flux is
computed by a bulk formula:

$$F_{sfc} = c_p \rho_1 C_H |u_1| (\theta_{sfc} - \theta_{atm1})$$
where $\theta_{sfc}$ is the potential temperature at the land or ocean surface,
$\theta_{atm1}$ is the potential temperature at the lowest atmospheric level,
$c_p$ is the specific heat, $C_H = 0.0015$ is the bulk transfer coefficient for
sensible heat, and $|u_1|$ is the near-surface atmospheric wind speed.
We assume that the potential temperature is defined with respect to the surface pressure,
so that $\theta_{sfc} = T_{sfc}$.
=#
coupler = CouplerState(Δt_coupled)
coupler_add_field!(coupler, :T_sfc_ocean, ocean.integrator.u.T_sfc; write_sim = ocean)
coupler_add_field!(coupler, :T_sfc_land, land.integrator.u.T_sfc; write_sim = land)
coupler_add_field!(coupler, :F_sfc, atmos.integrator.u.F_sfc; write_sim = atmos)
for (name, map) in pairs(maps)
    coupler_add_map!(coupler, name, map)
end

sim = AOLCoupledSim(atmos, ocean, land, coupler)


#=
## Coupled Time Integration

Finally, the execution sequence of the component models must be specified.
This is currently done explicitly with a combination of `step!`, `coupler_pull!`,
and `coupler_push!` methods. The `coupler_pull!` and `coupler_push!` methods
receive and send coupled field info from the coupler, respectively. They must
be written for each component simulation, and are simply collections of
`coupler_get` and `coupler_put!` methods for each component.

Here, the atmosphere steps forward first and then sends updated fields to
the coupler. The ocean and land (which are not coupled to each other) then
retreive the updated coupled information, advance and send their own updates
to the coupler.

Because the models exchange fluxes only at the coupled timestep, the surface flux
is accumulated over the coupled time-step
coupling time step,
`Δt_cpl`

$$F_{integ} = \int_{\Delta t_{coupler}} F_{sfc}  dt$$
where  $F_{integ}$ has units of $J m^{-2}$.
=#
function cpl_run(simulation::AOLCoupledSim)
    @info "Run model"
    (; atmos, ocean, land, coupler) = simulation
    Δt_coupled = coupler.Δt_coupled
    ## coupler stepping
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

## Run simulation
cpl_run(sim)

# ### References
# - [Antonelli & Rotunno 2007](https://journals.ametsoc.org/view/journals/atsc/64/12/2007jas2261.1.xml?tab_body=pdf)

## Post-processing
using JLD2 #hide
import Plots, ClimaCorePlots #hide

sol = sim.atmos.integrator.sol #hide
path = joinpath(@__DIR__, "output") #hide
mkpath(path) #hide
save(joinpath(path, "last_sim.jld2"), "coupled_sim", sim) #hide

Plots.GRBackend() #hide

anim = Plots.@animate for u in sol.u #hide
    Plots.contourf(u.Yc.ρθ ./ u.Yc.ρ) #hide
end #hide
Plots.mp4(anim, joinpath(path, "theta.mp4"), fps = 20) #hide

If2c = Operators.InterpolateF2C() #hide
anim = Plots.@animate for u in sol.u #hide
    Plots.contourf(If2c.(u.ρw) ./ u.Yc.ρ) #hide
end #hide

Plots.mp4(anim, joinpath(path, "vel_w.mp4"), fps = 20) #hide
anim = Plots.@animate for u in sol.u #hide
    Plots.contourf(u.Yc.ρuₕ ./ u.Yc.ρ) #hide
end #hide
Plots.mp4(anim, joinpath(path, "vel_u.mp4"), fps = 20) #hide
