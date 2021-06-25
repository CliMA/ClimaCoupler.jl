using ClimateMachine.ODESolvers
using ClimateMachine.ODESolvers: AbstractODESolver

export CplSolver

"""
    CplSolver(; component_list, coupler::CplState, coupling_dt, t0)

A time stepping like object for advancing a coupled system made up of a pre-defined
set of named components specified in `component_list`. Each component is a
balance law, discretization and timestepper collection. The coupler will
step them forward by a nsteps substeps to advance the coupled system
by a simulated time `coupling_dt`.

Components interact with the coupler during `pre_step()` and `post_step()` functions.
During a `pre_step()`, a component may get fields from the coupler name space.
A component may put fields into the coupler name space during the `post_step()` for
later use by other components. The `CplSolver` abstraction controls
 1. the outer time stepping sequencing of components
 2. the execution of actions mapping exports from one or more components to imports of
    other components through an intermediary coupler name space.
"""
mutable struct CplSolver{CL, FT, FL} <: AbstractODESolver
    "Named list of pre-defined components"
    component_list::CL
    "Coupler State"
    coupler::CplState
    "Coupling timestep"
    dt::FT
    "Start time - initializes or tries to restart"
    t0::FT
    "Current time"
    t::FT
    "elapsed number of steps"
    steps::Int
    "Flux Log"
    fluxlog::FL
end

function CplSolver(;
    component_list = component_list,
    coupler::CplState = coupler,
    coupling_dt = coupling_dt,
    t0 = t0,
    fluxlog = nothing,
)
    return CplSolver(component_list, coupler, coupling_dt, t0, t0, 0, fluxlog)
end

function ODESolvers.dostep!(Qtop, csolver::CplSolver, param, time::Real)

    for cpl_component in csolver.component_list

        # pre_step fetching imports goes here
        cpl_component.pre_step(csolver)
        component = cpl_component[:component_model]
        solve!(
            component.state,
            component.odesolver;
            numberofsteps = component.nsteps,
        )
        # post step pushing exports goes here
        cpl_component.post_step(csolver)

    end
    return nothing
end
