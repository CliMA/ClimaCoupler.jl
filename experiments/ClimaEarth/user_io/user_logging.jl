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
    checkpoint_callback(checkpointing_period::Dates.Period; active = true, start_date)

Create a schedule for checkpointing based on a time period.

This function returns an `EveryCalendarDtSchedule` object if `active` is true,
which can be used to schedule checkpoints at regular intervals defined by
`checkpointing_period`. If `active` is false, it returns `nothing`, effectively
disabling checkpointing.

# Arguments
- `checkpointing_period::Dates.Period`: The time period between checkpoints.
   This should be a `Dates.Period` object, e.g., `Dates.Month(1)` for monthly
   checkpoints.
- `active::Bool`: Whether checkpointing is active. Defaults to `true`.
- `start_date`: The start date of the simulation.

# Note

This function does not really return a callback, it returns a schedule that can
be used to call [`checkpoint_sims`](@ref) with [`maybe_checkpoint`](@ref).

# Example
```julia
callback = checkpoint_callback(Dates.Month(1); start_date = Date(2024, 1, 1))
```
"""
function checkpoint_callback(checkpointing_period::Dates.Period;
                             active = true,
                             start_date)
    active || return nothing
    return EveryCalendarDtSchedule(checkpointing_period; start_date)
end

"""
    maybe_checkpoint(cs, t, schedule)

Save a checkpoint, if it is time (as determined by `schedule`)

This function evaluates the given `schedule` at the current time `t` and
performs a checkpoint of the simulation state `cs` if the schedule dictates it.

# Arguments
- `cs`: The coupled simulation.
- `t`: The current time in the simulation.
- `schedule`: A callable object representing the checkpointing schedule.
  Typically produced with `checkpoint_callback`
"""
function maybe_checkpoint(cs, t, schedule)
    # EveryCalendarDtSchedule needs an integrator object and uses t from the object,
    # so we can just pass a namedtuple with t.
    schedule((; t)) && checkpoint_sims(cs, t)
    return nothing
end

maybe_checkpoint(cs, t, cb::Nothing) = nothing
