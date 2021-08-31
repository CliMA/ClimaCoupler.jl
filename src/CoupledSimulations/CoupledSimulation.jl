export CoupledSimulation, step!, run!

"""
    CoupledSimulation

An abstract type representing a coupled simulation.
"""
abstract type CouplerSimulation end

"""
    run!(::CoupledSimulation)

A simple outer timestepping loop for coupled system runs.

This will be formalized when the run! functionality for component
models is implemented so to have a consistent interface.
"""
function run!(sim::CoupledSimulation) end

"""
    step!(::CoupledSimulation, dt)

Advances a simulation by the coupling period `dt`.

This can be custom-written for different implementations of the
`CoupledSimulation` type. 
"""
function step!(sim::CoupledSimulation, dt) end