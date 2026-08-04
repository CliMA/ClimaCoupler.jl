#=
    Unit tests for ClimaCoupler TimeManager module
=#
import Test: @testset, @test, @test_logs, @test_throws
import Dates
import ClimaDiagnostics as CD
import ClimaUtilities.TimeManager: ITime
import ClimaCoupler: TimeManager

@testset "time_to_period" begin
    @test TimeManager.time_to_period("2months") == Dates.Month(2)
    @test TimeManager.time_to_period("10secs") == Dates.Millisecond(10_000)
    @test TimeManager.time_to_period("2.5hours") == Dates.Millisecond(9_000_000)
    @test_throws ErrorException TimeManager.time_to_period("never")
end

# A minimal integrator-like object, as built by `maybe_trigger_callback`
fake_integrator(t, Δt_cpl = 0.5) = (; t = Ref(t), tspan = (0.0, 100.0), Δt_cpl)

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

    TimeManager.maybe_trigger_callback(cb, fake_integrator(0.5))
    @test n_triggered[] == 0
    TimeManager.maybe_trigger_callback(cb, fake_integrator(2.0))
    @test n_triggered[] == 1

    never_cb = TimeManager.Callback(TimeManager.NeverSchedule(), cs -> n_triggered[] += 1)
    TimeManager.maybe_trigger_callback(never_cb, fake_integrator(2.0))
    @test n_triggered[] == 1
end

@testset "coupler_step_number" begin
    # Floating-point times
    @test TimeManager.coupler_step_number((;
        t = Ref(0.0),
        tspan = (0.0, 10.0),
        Δt_cpl = 2.0,
    )) == 0
    @test TimeManager.coupler_step_number((;
        t = Ref(2.0),
        tspan = (0.0, 10.0),
        Δt_cpl = 2.0,
    )) == 1
    @test TimeManager.coupler_step_number((;
        t = Ref(8.0),
        tspan = (0.0, 10.0),
        Δt_cpl = 2.0,
    )) == 4
    # Counted from the start of the current run, not of the simulation
    @test TimeManager.coupler_step_number((;
        t = Ref(8.0),
        tspan = (6.0, 10.0),
        Δt_cpl = 2.0,
    )) == 1
    # Robust to floating-point error (0.1 is not exactly representable)
    @test TimeManager.coupler_step_number((;
        t = Ref(sum(fill(0.1, 7))),
        tspan = (0.0, 10.0),
        Δt_cpl = 0.1,
    )) == 7

    # ITimes: dividing two ITimes gives an exact Rational
    itime(s) = ITime(s, epoch = Dates.DateTime(2010))
    @test TimeManager.coupler_step_number((;
        t = Ref(itime(0)),
        tspan = (itime(0), itime(100)),
        Δt_cpl = itime(20),
    )) == 0
    @test TimeManager.coupler_step_number((;
        t = Ref(itime(80)),
        tspan = (itime(0), itime(100)),
        Δt_cpl = itime(20),
    )) == 4
    @test TimeManager.coupler_step_number((;
        t = Ref(itime(80)),
        tspan = (itime(60), itime(100)),
        Δt_cpl = itime(20),
    )) == 1
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
        TimeManager.maybe_trigger_callback(cb, fake_integrator(0.5 * step))
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
        TimeManager.maybe_trigger_callback(cb, fake_integrator(1.0 * step, 1.0))
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

@testset "WalltimeReporter" begin
    @test TimeManager.compact_time_str(0.0) == "0 s"
    @test TimeManager.compact_time_str(59580.0) == "16 h 33 m"
    @test TimeManager.compact_time_str(2 * 86400.0) == "2 d"

    reporter = TimeManager.WalltimeReporter()
    fake_cs(t) = (;
        t = Ref(t),
        Δt_cpl = 400.0,
        tspan = (0.0, 86400.0),
        start_date = Dates.DateTime(2010),
    )

    # The first call reports progress but discards the wall time (compilation)
    @test_logs (:info, r"^Progress\n  time = 2010-01-01T00:06:40") reporter(fake_cs(400.0))
    @test reporter.wall_time_elapsed[] == 0.0

    # The second call reports timing estimates; its measurement is scaled by
    # (t - t_start) / (t - t_previous) = 2 to cover the pre-compilation steps
    sleep(0.1)
    @test_logs (:info, r"walltime remaining ≈ .*\n  sypd ≈ ") reporter(fake_cs(800.0))
    @test reporter.wall_time_elapsed[] >= 2 * 0.09

    # Later calls accumulate wall time without scaling
    elapsed_before = reporter.wall_time_elapsed[]
    sleep(0.1)
    @test_logs (:info, r"^Progress") reporter(fake_cs(1200.0))
    @test 0.09 <= reporter.wall_time_elapsed[] - elapsed_before <= 1.0
end
