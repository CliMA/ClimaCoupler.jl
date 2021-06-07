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

Components are registered with pre_step() and post_step() functions. The pre_step()
functions get fields for use in the component from the coupler name space.
The post_step() functions put fields for use by other components into
te coupler name space. 
The CplSolver abstraction controls
 1. the outer time stepping sequencing of components
 2. the excution of actions mapping exports from one or more components to imports of
    other components through an intermediary coupler name space.

Some notes -

For now components need to include slightly wasteful "shadow" variables for
accumulating boundary flux terms they compute across RK stages and across
timesteps. These are defined within the component balance law.
The shadown variable is a full 3d array because of the way the current 
infrastructure works. This can be tidied up later once design is settled.
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
