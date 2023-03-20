export CoupledSimulation, step!, run!, name

"""
    AbstractSimulation

An abstract type representing a model simulation.
"""
abstract type AbstractSimulation end

abstract type AbstractAtmosSimulation <: AbstractSimulation end
name(::AbstractAtmosSimulation) = :atmos

abstract type AbstractSurfaceSimulation <: AbstractSimulation end

abstract type AbstractOceanSimulation <: AbstractSurfaceSimulation end
name(::AbstractOceanSimulation) = :ocean

abstract type AbstractLandSimulation <: AbstractSurfaceSimulation end
name(::AbstractLandSimulation) = :land

abstract type AbstractCoupledSimulation <: AbstractSimulation end
name(::AbstractCoupledSimulation) = :coupled

"""
    CoupledSimulation

A subtype of the abstract type `AbstractCoupledSimulation` representing a model simulation.
"""
# struct CoupledSimulation{CS, S, CPL, L, C} <: AbstractCoupledSimulation
#     "The coupled time-stepping scheme"
#     coupler_solver::CS
#     "The component simulations"
#     simulations::S
#     "The coupler"
#     coupler::CPL
#     "Diagnostic logger"
#     logger::L
#     "Clock"
#     clock::C
# end


mutable struct CouplerState{FT, CF, RO}
    # A dictionary of fields added to the coupler
    coupled_fields::CF
    # A dictionary of remap operators between components
    remap_operators::RO
    # The coupled timestep size
    Δt_coupled::FT
end


struct CoupledSimulation{CS, S, CPL, L, C} <: AbstractCoupledSimulation
    "Communication (MPI/GPU) context"
    comms_context
    "Calendar dates"
    dates

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

mutable struct cs
    comms_context
    parsed_args
    component_sims
    exchange_masks # (on coupler grid)
    exchange_fields
    clock
    calendar
    diagnostics # callback
    conservation_checks # callback
    flux_calculator
    regrid_maps
    metadata # boundary_space? or func

    # callbacks # TBD
    # running
    # initialized
end

mutable struct clock
    timestepper
    tspan
    dt
    current_time
end

mutable struct calendar
    start_date
    current_date
    first_day_of_month
end

mutable struct component_sim
    model
    integrator
    parameters
    mask
end


"""
    run!(::CoupledSimulation)

A simple outer timestepping loop for coupled system runs.

This will be formalized when the `run!` functionality for component
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

Advances a simulation `sim` by `dt`.

Note that `dt` is not necessarily the simulation's timestep length;
a simuation could take several shorter steps that total to `dt`.
"""
function step!(sim::AbstractSimulation, dt) end


# each component model should define this
get_temperature(::AbstractSimulation)= nothing
get_global_energy(::AbstractSimulation)= nothing
coupler_get(::AbstractSimulation) = nothing #?
coupler_put(::AbstractSimulation)= nothing #?


# each experiment should define this
# CplFieldInfo(name) # for all exchsnge fields
# timestepping_order() # to be replaced by Timestepper module when concurrent coupling implemented