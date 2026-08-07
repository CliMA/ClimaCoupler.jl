#
#   Unit tests for ClimaCoupler TimeManager module
#
import Test: @testset, @test, @test_logs, @test_throws
import Dates
import ClimaDiagnostics as CD
import ClimaCoupler: TimeManager

@testset "time_to_period" begin
    @test TimeManager.time_to_period("2months") == Dates.Month(2)
    @test TimeManager.time_to_period("10secs") == Dates.Millisecond(10_000)
    @test TimeManager.time_to_period("2.5hours") == Dates.Millisecond(9_000_000)
    @test_throws ErrorException TimeManager.time_to_period("never")
end

# A minimal coupled simulation, of which `maybe_trigger_callback` only reads `t` and `step`
fake_cs(t, step = 1) = (; t = Ref(t), step = Ref(step))

# A schedule that is never true, but records how many times it was asked
struct CountingSchedule <: CD.Schedules.AbstractSchedule
    n_calls::Base.RefValue{Int}
end
(schedule::CountingSchedule)(_) = (schedule.n_calls[] += 1; false)
CD.Schedules.short_name(::CountingSchedule) = "counting"
CD.Schedules.long_name(::CountingSchedule) = "counting"

@testset "Callback triggering" begin
    n_triggered = Ref(0)
    cb = TimeManager.Callback(integrator -> integrator.t > 1.0, cs -> n_triggered[] += 1)

    TimeManager.maybe_trigger_callback(cb, fake_cs(0.5))
    @test n_triggered[] == 0
    TimeManager.maybe_trigger_callback(cb, fake_cs(2.0))
    @test n_triggered[] == 1

    never_cb = TimeManager.Callback(TimeManager.NeverSchedule(), cs -> n_triggered[] += 1)
    TimeManager.maybe_trigger_callback(never_cb, fake_cs(2.0))
    @test n_triggered[] == 1
end

@testset "NeverSchedule" begin
    schedule = TimeManager.NeverSchedule()
    @test schedule((; t = 0.0, step = 1)) == false
    @test schedule((; t = 100.0, step = 12)) == false
    @test CD.Schedules.short_name(schedule) == "never"
    @test CD.Schedules.long_name(schedule) == "never"
end

@testset "OnceSchedule" begin
    schedule = TimeManager.OnceSchedule()
    @test schedule((; t = 1.0, step = 1)) == true
    @test schedule((; t = 0.0, step = 0)) == false
    @test schedule((; t = 2.0, step = 2)) == false
    @test CD.Schedules.short_name(schedule) == "once"
end

@testset "PowerOfTwoSchedule" begin
    schedule = TimeManager.PowerOfTwoSchedule()
    triggered = filter(step -> schedule((; t = float(step), step)), 0:20)
    @test triggered == [1, 2, 4, 8, 16]

    @test CD.Schedules.short_name(schedule) == "pow2"
    @test CD.Schedules.long_name(schedule) == "every power-of-two iteration"

    # Callbacks see a step number, so the schedule works end-to-end
    n_triggered = Ref(0)
    cb = TimeManager.Callback(schedule, cs -> n_triggered[] += 1)
    for step in 1:8
        TimeManager.maybe_trigger_callback(cb, fake_cs(0.5 * step, step))
    end
    @test n_triggered[] == 4 # steps 1, 2, 4, 8
end

@testset "OrSchedule" begin
    never = TimeManager.NeverSchedule()
    always = CD.Schedules.EveryStepSchedule()
    integrator = (; t = 3.0, step = 3)

    @test TimeManager.OrSchedule(never, never)(integrator) == false
    @test TimeManager.OrSchedule(never, always)(integrator) == true
    @test TimeManager.OrSchedule(always, never)(integrator) == true
    @test TimeManager.OrSchedule(never, never, always)(integrator) == true
    @test TimeManager.OrSchedule(never)(integrator) == false

    @test_throws ErrorException TimeManager.OrSchedule()

    # Every schedule must be called, even once the answer is known: stateful
    # schedules only advance when they return true, so skipping one would make
    # it fire late. Guards against short-circuiting `any`.
    counter = CountingSchedule(Ref(0))
    # `always` returns true first, but `counter` must still be queried
    @test TimeManager.OrSchedule(always, counter)(integrator) == true
    @test counter.n_calls[] == 1

    # The concrete regression: a stateful schedule that is skipped fires late
    n_triggered = Ref(0)
    or_schedule = TimeManager.OrSchedule(
        TimeManager.PowerOfTwoSchedule(),
        CD.Schedules.EveryDtSchedule(4.0),
    )
    cb = TimeManager.Callback(or_schedule, cs -> n_triggered[] += 1)
    triggered_steps = Int[]
    for step in 1:8
        n_before = n_triggered[]
        TimeManager.maybe_trigger_callback(cb, fake_cs(1.0 * step, step))
        n_triggered[] > n_before && push!(triggered_steps, step)
    end
    # steps 1, 2, 4, 8 from PowerOfTwoSchedule; 4 and 8 also from EveryDtSchedule.
    # Step 5 must NOT trigger: if the EveryDtSchedule had been skipped at step 4,
    # it would fire there instead.
    @test triggered_steps == [1, 2, 4, 8]

    # Names are composed from the children
    named = TimeManager.OrSchedule(TimeManager.PowerOfTwoSchedule(), never)
    @test CD.Schedules.short_name(named) == "pow2_or_never"
    @test CD.Schedules.long_name(named) == "every power-of-two iteration or never"
end

@testset "calendar_dt_schedule" begin
    start_date = Dates.DateTime(2010)

    # A fresh simulation is unchanged: nothing fires at t = 0, first firing at day 1
    fresh = TimeManager.calendar_dt_schedule("1days", start_date, 0.0)
    @test !fresh((; t = 0.0))
    @test fresh((; t = 86400.0))

    # A simulation restarted at day 2 stays on the calendar boundaries of the
    # original run: nothing fires one coupling step after the restart, the first
    # firing is at day 3 exactly, and not again one step later
    restarted = TimeManager.calendar_dt_schedule("1days", start_date, 172800.0)
    @test !restarted((; t = 172980.0))
    @test restarted((; t = 259200.0))
    @test !restarted((; t = 259380.0))

    # A Dates.Period is accepted as well as a time string
    periodic = TimeManager.calendar_dt_schedule(Dates.Day(1), start_date, 0.0)
    @test periodic((; t = 86400.0))
end

@testset "walltime_schedule" begin
    start_date = Dates.DateTime(2010)
    periodic = CD.Schedules.EveryCalendarDtSchedule(
        TimeManager.time_to_period("1days");
        start_date,
    )

    # Without walltime_debug, reporting is purely periodic
    @test TimeManager.walltime_schedule("1days", false, start_date) == periodic
    @test isnothing(TimeManager.walltime_schedule("never", false, start_date))

    # With walltime_debug, the powers of two are added on top
    schedule = TimeManager.walltime_schedule("1days", true, start_date)
    @test schedule isa TimeManager.OrSchedule
    # `time_to_period` builds the period in milliseconds, hence the name
    @test CD.Schedules.short_name(schedule) == "86400000millis_or_pow2"
    # ... reports on step 2, which is not a multiple of a day
    @test schedule((; t = 800.0, step = 2)) == true
    # ... and still reports periodically, on a step that is not a power of two
    @test schedule((; t = 3 * 86400.0, step = 3)) == true

    # If walltime_dt is "never", the powers of two are the only reporting
    schedule = TimeManager.walltime_schedule("never", true, start_date)
    @test schedule == TimeManager.PowerOfTwoSchedule()
    @test schedule((; t = 800.0, step = 2)) == true
    @test schedule((; t = 1200.0, step = 3)) == false
end

@testset "WalltimeReporter" begin
    @test TimeManager.compact_time_str(0.0) == "0 s"
    @test TimeManager.compact_time_str(59580.0) == "16 h 33 m"
    @test TimeManager.compact_time_str(2 * 86400.0) == "2 d"

    reporter = TimeManager.WalltimeReporter()
    reporter_cs(t) = (;
        t = Ref(t),
        step = Ref(round(Int, t / 400.0)),
        Δt_cpl = 400.0,
        tspan = (0.0, 86400.0),
        start_date = Dates.DateTime(2010),
    )

    # The first call reports progress but discards the wall time (compilation)
    @test_logs (:info, r"^Progress\n  time = 2010-01-01T00:06:40") reporter(
        reporter_cs(400.0),
    )
    @test reporter.wall_time_elapsed[] == 0.0

    # The second call reports timing estimates; its measurement is scaled by
    # (t - t_start) / (t - t_previous) = 2 to cover the pre-compilation steps
    sleep(0.1)
    @test_logs (:info, r"walltime remaining ≈ .*\n  sypd ≈ ") reporter(reporter_cs(800.0))
    @test reporter.wall_time_elapsed[] >= 2 * 0.09

    # Later calls accumulate wall time without scaling
    elapsed_before = reporter.wall_time_elapsed[]
    sleep(0.1)
    @test_logs (:info, r"^Progress") reporter(reporter_cs(1200.0))
    @test 0.09 <= reporter.wall_time_elapsed[] - elapsed_before <= 1.0
end
