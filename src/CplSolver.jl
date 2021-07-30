using ClimateMachine.ODESolvers
using ClimateMachine.ODESolvers: AbstractODESolver
using Base.Threads

export SerialCplSolver, ConcurrentCplSolver

"""
    AbstractCplSolver

An abstract type for coupled solver schemes.

At present, this extends ClimateMachine's AbstractODESolver interface.
Concrete implementations currently include `SerialCplSolver` and `ConcurrentCplSolver`. 
Coupling schemes such as implicit schemes or higher order explicit schemes may need
to store additional info such as convergence tolerances or substage history.
"""
# TODO: may not tie this to the AbstractODESolver in the future. Depends on outer time loop mechanism
abstract type AbstractCplSolver <: AbstractODESolver end

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

"""
    SerialCplSolver(; component_list, coupler, coupling_dt, t0, fluxlog)

A sequential-leapfrog coupling scheme.
"""
mutable struct SerialCplSolver{CL, FT, FL} <: AbstractCplSolver
    component_list::CL # Named list of pre-defined components
    coupler::CplState # Coupler State
    dt::FT # Coupling timestep
    t0::FT # Start time - initializes or tries to restart
    t::FT # Current time
    steps::Int # elapsed number of steps
    fluxlog::FL # Flux Log
end

function SerialCplSolver(;
    component_list = component_list,
    coupler::CplState = coupler,
    coupling_dt = coupling_dt,
    t0 = t0,
    fluxlog = nothing,
)
    return SerialCplSolver(component_list, coupler, coupling_dt, t0, t0, 0, fluxlog)
end

function ODESolvers.dostep!(Qtop, csolver::SerialCplSolver, param, time::Real)
    for cpl_component in csolver.component_list

        # pre_step fetching imports goes here
        cpl_component.pre_step(csolver)
        component = cpl_component[:component_model]
        solve!(
            component.state,
            component.odesolver,
            param;
            numberofsteps = component.nsteps,
            callbacks = component.callbacks,
        )
        # post step pushing exports goes here
        cpl_component.post_step(csolver)

    end
    return nothing
end

"""
    ConcurrentCplSolver(; component_list, coupler, coupling_dt, t0, fluxlog)

A parallel explicit coupling scheme using Julia's multithreading.
"""
mutable struct ConcurrentCplSolver{CL, FT, FL} <: AbstractCplSolver
    component_list::CL # Named list of pre-defined components
    coupler::CplState # Coupler State
    dt::FT # Coupling timestep
    t0::FT # Start time - initializes or tries to restart
    t::FT # Current time
    steps::Int # elapsed number of steps
    fluxlog::FL # Flux Log
end

function ConcurrentCplSolver(;
    component_list = component_list,
    coupler::CplState = coupler,
    coupling_dt = coupling_dt,
    t0 = t0,
    fluxlog = nothing,
)
    return ConcurrentCplSolver(component_list, coupler, coupling_dt, t0, t0, 0, fluxlog)
end

function ODESolvers.dostep!(Qtop, csolver::ConcurrentCplSolver, param, time::Real)
    @threads for cpl_component in csolver.component_list
        # pre_step fetching imports goes here
        cpl_component.pre_step(csolver)
    end
    @threads for cpl_component in csolver.component_list
            component = cpl_component[:component_model]
            solve!(
            component.state,
            component.odesolver,
            param;
            numberofsteps = component.nsteps,
            callbacks = component.callbacks,
        )
    end
    @threads for cpl_component in csolver.component_list
        # post step pushing exports goes here
        cpl_component.post_step(csolver)
    end
    return nothing
end
