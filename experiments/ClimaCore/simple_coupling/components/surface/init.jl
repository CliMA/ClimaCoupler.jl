struct ThermalSlab{FT} <: AbstractComponentModel end

# Slab RHS
"""
    get_rhs!(::ThermalSlab{FT})

Slab layer equation
    lnd d(T_sfc)/dt = - F_accumulated / h_lnd
    where
        F_accumulated = F_integrated / Δt_coupler
"""
function get_rhs!(::ThermalSlab{FT})
    function rhs!(du, u, cache, t)
        @. du.T = (-cache.flux) / cache.params.h
    end
end

# column init
function init(tag::ThermalSlab{FT}; timestepping = (1,1,1), space = nothing) where {FT} #; tspan, dt, saveat, space, ocean_mask, stepper = Euler()) where {FT}

    tspan, dt, saveat = timestepping

    # column parameters
    params = (
        h = FT(0.5), # depth of slab layer [m]
        T_0 = FT(260.0), # initial condition o
        )

    # initialize prognostic variables, either as ClimaCore's Field objects or as Arrays
    Y_0 = Fields.FieldVector(; T = Fields.ones(space) .* params.T_0,)

    cache = (
        params = params,
        flux = Fields.zeros(space),
    ) # cache needed for ODEProblem

    problem = OrdinaryDiffEq.ODEProblem(get_rhs!(tag), Y_0, tspan, cache)

    stepper = Euler()
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    ComponentSimulation(tag, Y_init = Y_0, integrator = integrator)
end