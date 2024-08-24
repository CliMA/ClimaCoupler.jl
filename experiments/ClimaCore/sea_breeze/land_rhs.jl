# # Land Model

import DiffEqCallbacks
import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS

# Load coupled simulation code
include("../CoupledSims/coupled_sim.jl")

#=
## Slab Land ODE
For our land component, we solve a simple slab land ODE:

$$\rho_l c_l H_l \partial_t T_{lnd} = - F_{integ} / \Delta t_{coupler}$$
- where $\rho_l = 1500$ kg m $^{-3}$, $c_l=800$ J K $^{-1}$ kg $^{-1}$, $H_l=1$ m are the density, specific heat and depth of the land slab,
- and $F_{integ}$ is the integrated surface fluxes in time.
=#

# ## Model Code
function lnd_rhs!(du, u, (parameters, F_accumulated), t)
    """
    Slab layer equation
        d(T_lnd)/dt = - (F_accumulated + G) / (h_lnd * ρ_lnd * c_lnd)
        where
            F_accumulated = F_integrated / Δt_coupler
    """
    (; lnd_h, lnd_ρ, lnd_c) = parameters
    (; T_sfc) = du

    @. T_sfc = (-F_accumulated) / (lnd_h * lnd_ρ * lnd_c)
end

## set up domain
function hspace_1D(xlim = (-π, π), npoly = 0, helem = 10)
    FT = Float64

    domain =
        CC.Domains.IntervalDomain(CC.Geometry.XPoint{FT}(xlim[1]), CC.Geometry.XPoint{FT}(xlim[2]), periodic = true)
    mesh = CC.Meshes.IntervalMesh(domain; nelems = helem)
    topology = CC.Topologies.IntervalTopology(mesh)

    ## Finite Volume Approximation: Gauss-Lobatto with 1pt per element
    quad = CC.Spaces.Quadratures.GL{npoly + 1}()
    space = CC.Spaces.SpectralElementSpace1D(topology, quad)

    return space
end

## init simulation
function lnd_init(; xmin = -1000, xmax = 1000, helem = 20, npoly = 0)

    ## construct domain spaces - get only surface layer (NB: z should be zero, not z = first central height)
    space = hspace_1D((xmin, xmax), npoly, helem)
    coords = CC.Fields.coordinate_field(space)
    domain = space

    ## initial condition
    T_sfc = map(coords) do coord
        T_sfc = 283.0
    end

    ## prognostic variable
    Y = CC.Fields.FieldVector(T_sfc = T_sfc)

    return Y, domain
end

# ## Coupled Land Wrappers
## Land Simulation - later to live in ClimaLand
struct LandSim <: AbstractLandSim
    integrator::Any
end

function LandSim(Y_init, t_start, dt, t_end, timestepper, p, saveat, callbacks = DiffEqCallbacks.CallbackSet())
    ode_algo = CTS.ExplicitAlgorithm(timestepper)
    ode_function = CTS.ClimaODEFunction(T_exp! = lnd_rhs!)

    problem = SciMLBase.ODEProblem(ode_function, Y_init, (t_start, t_end), p)
    lnd_integ = SciMLBase.init(problem, ode_algo, dt = dt, saveat = saveat, adaptive = false, callback = callbacks)

    return LandSim(lnd_integ)
end

function coupler_push!(coupler::CouplerState, land::LandSim)
    coupler_put!(coupler, :T_sfc_land, land.integrator.u.T_sfc, land)
end

function coupler_pull!(land::LandSim, coupler::CouplerState)
    coupler_get!(land.integrator.p.F_sfc, coupler, :F_sfc, land)
    land.integrator.p.F_sfc ./= coupler.Δt_coupled
end
