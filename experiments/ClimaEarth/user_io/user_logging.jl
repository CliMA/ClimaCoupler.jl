import ClimaCoupler: Checkpointer, Interfacer

import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule

"""
    Base.show(io::IO, dict::Dict)

This prints the keys and values of a Dict in sorted order.
"""
function Base.show(io::IO, dict::Dict)
    for k in sort!(collect(keys(dict)))
        println(io, " $k => $(dict[k])")
    end
end

# user callbacks
"""
    checkpoint_sims(cs::CoupledSimulation)

This is a callback function that checkpoints all simulations defined in the current coupled simulation.

Checkpoints are saved in the `cs.dirs.checkpoints` folder.
"""
function checkpoint_sims(cs::Interfacer.CoupledSimulation, t)
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            Checkpointer.checkpoint_model_state(sim, cs.comms_ctx, Int(round(t)), output_dir = cs.dirs.checkpoints)
        end
    end
end

"""
    periodic_callback(period::Dates.Period; start_date)

Create a schedule for  based on a time period.

This function returns an `EveryCalendarDtSchedule` object which can be used to
schedule checkpoints at regular intervals defined by `period`..

# Arguments
- `period::Dates.Period`: The time period between calls.
   This should be a `Dates.Period` object, e.g., `Dates.Month(1)` for monthly
   calls.
- `start_date`: The start date of the simulation.

# Note

This function does not really return a callback, it returns a schedule that can
be used to call the desired function

# Example
```julia
callback = periodic_callback(Dates.Month(1); start_date = Date(2024, 1, 1))
```
"""
function periodic_callback(period::Dates.Period; start_date)
    return EveryCalendarDtSchedule(period; start_date)
end
