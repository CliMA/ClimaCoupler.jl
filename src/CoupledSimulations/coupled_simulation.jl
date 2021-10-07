export CoupledSimulation, step!, run!

"""
    CoupledSimulation

An abstract type representing a coupled simulation.
"""
abstract type CoupledSimulation end

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
function step!(sim, dt) end