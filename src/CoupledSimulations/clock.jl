"""
    Clock{T}

Manages a simulation's time information.
"""
mutable struct Clock{T}
    time::T         # current simulation time
    dt::T           # simulation timestep
    stop_time::T    # simulation end time
end

tick!(clock::Clock) = (clock.time += clock.dt)

stop_time_exceeded(clock::Clock) = (clock.time >= clock.stop_time)
