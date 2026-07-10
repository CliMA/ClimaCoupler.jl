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
        isnothing(months) &&
            error("$(s) has to be of the form <NUM>months, e.g. 2months for 2 months")
        return Dates.Month(parse(Int, first(months)))
    else
        # Milliseconds to support fractional seconds
        return Dates.Millisecond(1000 * time_to_seconds(s))
    end
end

"""
    simulated_years_per_day(t_start, t_end, walltime)

Compute the simulated years per walltime day for a simulation spanning from
`t_start` to `t_end` (in seconds), assuming that the simulation took `walltime`
(in seconds) to run.
"""
simulated_years_per_day(t_start, t_end, walltime) =
    float(t_end - t_start) / walltime / 365.25

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

Check if it is time to call `callback`, if yes, call its function on `cs`.
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

"""
    WalltimeReporter()

A callable object that logs the progress of a coupled simulation.

A `WalltimeReporter` is meant to be used as the function of a `Callback`: it is
called with the `CoupledSimulation` as its only argument and reads the current
time, start date, time span, and coupling time step from it.

The wall time of the first call is dominated by compilation, so it is not
measured; instead, the first measurement (between the first and second calls) is 
scaled up to estimate the compilation-free wall time of the steps before it.
"""
struct WalltimeReporter
    n_calls::Base.RefValue{Int}
    wall_time_last::Base.RefValue{Float64}
    sim_time_last::Base.RefValue{Float64}
    wall_time_elapsed::Base.RefValue{Float64}
end
WalltimeReporter() = WalltimeReporter(Ref(0), Ref(0.0), Ref(0.0), Ref(0.0))

function (reporter::WalltimeReporter)(cs)
    # t, t_start, t_end, Δt_cpl, and date from the coupled simulation
    t = float(cs.t[])
    t_start, t_end = float.(cs.tspan)
    Δt_cpl = float(cs.Δt_cpl)
    date = cs.start_date + Dates.Second(round(Int, t))

    # compute the number of steps completed, total number of steps, and percent complete
    n_steps = max(1, round(Int, (t - t_start) / Δt_cpl))
    n_steps_total = ceil(Int, (t_end - t_start) / Δt_cpl) #TODO: This doesn't support t_end = Inf
    percent_complete = round(100 * (t - t_start) / (t_end - t_start); digits = 1)

    # begin info message
    msg = """
        Progress
          time = $date ($(compact_time_str(t - t_start)))
          step = $n_steps ($percent_complete%)"""

    reporter.n_calls[] += 1
    if reporter.n_calls[] > 1
        Δt_wall = time() - reporter.wall_time_last[]
        if reporter.wall_time_elapsed[] == 0
            # Scale the first measurement to also cover the steps before the
            # previous call, whose wall time was discarded as compilation
            Δt_wall *= (t - t_start) / (t - reporter.sim_time_last[])
        end
        reporter.wall_time_elapsed[] += Δt_wall

        wall_time_remaining =
            reporter.wall_time_elapsed[] / n_steps * (n_steps_total - n_steps)
        finish_date =
            trunc(Dates.now() + Dates.Second(round(Int, wall_time_remaining)), Dates.Second)
        sypd = simulated_years_per_day(t_start, t, reporter.wall_time_elapsed[])
        msg *= "\n  walltime remaining ≈ $(compact_time_str(wall_time_remaining)) ($finish_date)"
        msg *= "\n  sypd ≈ $(round(sypd; sigdigits = 4))"
    end
    reporter.sim_time_last[] = t
    reporter.wall_time_last[] = time()
    @info msg
    return nothing
end

"""
    compact_time_str(seconds)

Return a compact string for a duration given in `seconds` using the two largest nonzero units. 

Example: `compact_time_str(59580.0) == "16 h 33 m"`.
"""
function compact_time_str(seconds)
    remaining = round(Int, seconds)
    remaining < 1 && return "0 s"
    parts = String[]
    for (unit, unit_length) in
        (("y", 31557600), ("d", 86400), ("h", 3600), ("m", 60), ("s", 1))
        n, remaining = divrem(remaining, unit_length)
        n > 0 && push!(parts, "$n $unit")
        length(parts) == 2 && break
    end
    return join(parts, " ")
end

end
