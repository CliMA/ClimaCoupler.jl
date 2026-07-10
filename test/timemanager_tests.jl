#=
    Unit tests for ClimaCoupler TimeManager module
=#
import Test: @testset, @test, @test_logs, @test_throws
import Dates
import ClimaCoupler: TimeManager

@testset "time_to_period" begin
    @test TimeManager.time_to_period("2months") == Dates.Month(2)
    @test TimeManager.time_to_period("10secs") == Dates.Millisecond(10_000)
    @test TimeManager.time_to_period("2.5hours") == Dates.Millisecond(9_000_000)
    @test_throws ErrorException TimeManager.time_to_period("never")
end

@testset "Callback triggering" begin
    n_triggered = Ref(0)
    cb = TimeManager.Callback(integrator -> integrator.t > 1.0, cs -> n_triggered[] += 1)

    TimeManager.maybe_trigger_callback(cb, (; t = Ref(0.5)))
    @test n_triggered[] == 0
    TimeManager.maybe_trigger_callback(cb, (; t = Ref(2.0)))
    @test n_triggered[] == 1

    never_cb = TimeManager.Callback(TimeManager.NeverSchedule(), cs -> n_triggered[] += 1)
    TimeManager.maybe_trigger_callback(never_cb, (; t = Ref(2.0)))
    @test n_triggered[] == 1
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
