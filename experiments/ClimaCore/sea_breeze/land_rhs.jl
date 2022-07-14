# # Land Model

#=
## Slab Land ODE
For our land component, we solve a simple slab land ODE:

$$\rho_l c_l H_l \partial_t T_{lnd} = - F_{integ} / \Delta t_{coupler}$$
- where $\rho_l = 1500$ kg m $^{-3}$, $c_l=800$ J K $^{-1}$ kg $^{-1}$, $H_l=1$ m are the density, specific heat and depth of the land slab,
- and $F_{integ}$ is the integrated surface fluxes in time.
=#

function lnd_rhs!(du, u, (parameters, F_accumulated), t)
    """
    Slab layer equation
        d(T_lnd)/dt = - (F_accumulated + G) / (h_lnd * ρ_lnd * c_lnd)
        where 
            F_accumulated = F_integrated / Δt_coupler
    """
    @unpack lnd_h, lnd_ρ, lnd_c = parameters
    @unpack T_sfc = du

    @. T_sfc = (-F_accumulated) / (lnd_h * lnd_ρ * lnd_c)
end

## set up domain
function hspace_1D(xlim = (-π, π), npoly = 0, helem = 10)
    FT = Float64

    domain = Domains.IntervalDomain(Geometry.XPoint{FT}(xlim[1]) .. Geometry.XPoint{FT}(xlim[2]), periodic = true)
    mesh = Meshes.IntervalMesh(domain; nelems = helem)
    topology = Topologies.IntervalTopology(mesh)

    ## Finite Volume Approximation: Gauss-Lobatto with 1pt per element
    quad = Spaces.Quadratures.GL{npoly + 1}()
    space = Spaces.SpectralElementSpace1D(topology, quad)

    return space
end

## init simulation
function lnd_init(; xmin = -1000, xmax = 1000, helem = 20, npoly = 0)

    ## construct domain spaces - get only surface layer (NB: z should be zero, not z = first central height)
    space = hspace_1D((xmin, xmax), npoly, helem)
    coords = Fields.coordinate_field(space)
    domain = space

    ## initial condition
    T_sfc = map(coords) do coord
        T_sfc = 273.0
    end

    ## prognostic variable
    Y = Fields.FieldVector(T_sfc = T_sfc)

    return Y, domain
end

# ## Coupled Land Wrappers
## Land Simulation - later to live in ClimaLSM
struct LandSimulation <: ClimaCoupler.AbstractLandSimulation
    integrator::Any
end

function LandSimulation(Y_init, t_start, dt, t_end, timestepper, p, saveat, callbacks = CallbackSet())
    lnd_prob = ODEProblem(lnd_rhs!, Y_init, (t_start, t_end), p)
    lnd_integ = init(lnd_prob, timestepper, dt = dt, saveat = saveat, callback = callbacks)
    return LandSimulation(lnd_integ)
end

function ClimaCoupler.coupler_push!(coupler::ClimaCoupler.CouplerState, land::LandSimulation)
    coupler_put!(coupler, :T_sfc_land, land.integrator.u.T_sfc, land)
end

function ClimaCoupler.coupler_pull!(land::LandSimulation, coupler::ClimaCoupler.CouplerState)
    coupler_get!(land.integrator.p.F_sfc, coupler, :F_sfc, land)
    land.integrator.p.F_sfc ./= coupler.Δt_coupled
end
