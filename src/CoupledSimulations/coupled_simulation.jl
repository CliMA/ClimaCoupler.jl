export CoupledSimulation, step!, run!, name

"""
    AbstractSimulation

An abstract type representing a model simulation.
"""
abstract type AbstractSimulation end

abstract type AbstractAtmosSimulation <: AbstractSimulation end
name(::AbstractAtmosSimulation) = :atmos

abstract type AbstractOceanSimulation <: AbstractSimulation end
name(::AbstractOceanSimulation) = :ocean

abstract type AbstractLandSimulation <: AbstractSimulation end
name(::AbstractLandSimulation) = :land

abstract type AbstractCoupledSimulation <: AbstractSimulation end
name(::AbstractCoupledSimulation) = :coupled

struct CoupledSimulation{CS, S, CPL, L, C} <: AbstractCoupledSimulation
    "The coupled time-stepping scheme"
    coupler_solver::CS
    "The component simulations"
    simulations::S
    "The coupler"
    coupler::CPL
    "Diagnostic logger"
    logger::L
    "Clock"
    clock::C
end



"""
    run!(::CoupledSimulation)

A simple outer timestepping loop for coupled system runs.

This will be formalized when the run! functionality for component
models is implemented so to have a consistent interface.
"""
function run!(sim::CoupledSimulation)
    clock = sim.clock
    while !stop_time_exceeded(clock)
        step!(sim, clock.dt)
        tick!(clock)
    end
end

"""
    step!(sim, dt)

Advances a simulation by `dt`.

Note that `dt` is not necessarily the simulation's timestep length;
a simuation could take several shorter steps that total to `dt`.
"""
function step!(sim::AbstractSimulation, dt) end
