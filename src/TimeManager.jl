"""
    TimeManager

This module facilitates calendar functions and temporal interpolations
of data.
"""
module TimeManager

import Dates
import ..Interfacer
import ..Utilities: time_to_seconds

"""
    time_to_period(s::String)

Convert a string to a `Dates.Period` object.

The string has to have format `<NUM><unit>`, where `<NUM>` is a number (integer
or floating point) and `<unit>` is one of `secs`, `mins`, `hours`, `days`, or
`months`.

# Arguments
- `s::String`: The string to convert to a `Dates.Period`.

# Returns
- A `Dates.Period` object representing the time period specified in the string.

# Examples
```julia
julia> time_to_period("2months")
2 months

julia> time_to_period("10secs")
10000 milliseconds

julia> time_to_period("2.5hours")
9000000 milliseconds
```
"""
function time_to_period(s::String)
    if occursin("months", s)
        months = match(r"^(\d+)months$", s)
        isnothing(months) && error("$(s) has to be of the form <NUM>months, e.g. 2months for 2 months")
        return Dates.Month(parse(Int, first(months)))
    else
        # Milliseconds to support fractional seconds
        return Dates.Millisecond(1000 * time_to_seconds(s))
    end
end

"""
    Callback

A small struct containing a schedule, and the function to be executed if the
schedule is true.

A `schedule` is a callable object (ie, a function) that takes an integrator-type
of object and returns true or false.

The function `func` calls the coupled state `cs`.
"""
struct Callback{SCHEDULE, FUNC}
    schedule::SCHEDULE
    func::FUNC
end

"""
    maybe_trigger_callback(callback, cs)

Check if it time to call `callback`, if yes, call its function on `cs`.
"""
function maybe_trigger_callback(callback, cs)
    t = cs.t[]
    callback.schedule((; t)) && callback.func(cs)
    return nothing
end

"""
    callback!(cs)

Run the callbacks that need to be run.

This is marked as mutating function because in general callbacks are going to change
the state of the simulation
"""
function callbacks!(cs)
    foreach(c -> TimeManager.maybe_trigger_callback(c, cs), cs.callbacks)
    return nothing
end


"""
A schedule that is never true. Useful to disable something.

TODO: Move this to ClimaUtilities once we move Schedules there
"""
struct NeverSchedule end
function (::NeverSchedule)(args...)
    return false
end

end
