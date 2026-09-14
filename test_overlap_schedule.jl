# Exercise the overlapped slow-surface schedule without running a model.
#
# The launch/join logic assumes a window: launch at the end of one coupling step,
# join at the top of a later one. When dt_ocean == dt_cpl (k = 1) those collapse
# onto consecutive steps, which is legal and needs checking.
#
#   julia --project=experiments/CMIP test_overlap_schedule.jl

using ClimaCoupler
import Dates
import ClimaUtilities.TimeManager: ITime, date
const Interfacer = ClimaCoupler.Interfacer
const FieldExchanger = ClimaCoupler.FieldExchanger
const SimCoordinator = ClimaCoupler.SimCoordinator

mutable struct StubOcean <: Interfacer.AbstractOceanSimulation
    clock::Float64
    dt::Float64
    nsteps::Int
end
mutable struct StubIce <: Interfacer.AbstractSeaIceSimulation
    clock::Float64
    dt::Float64
    nsteps::Int
end
const Stub = Union{StubOcean, StubIce}

Interfacer.sim_dt(s::Stub) = s.dt
Interfacer.will_step(s::Stub, t::Float64) = (Float64(t) - s.clock) >= s.dt
# `t::Float64` is required: Interfacer defines step!(::AbstractComponentSimulation,
# ::Float64), so an untyped `t` here would be ambiguous rather than more specific.
function Interfacer.step!(s::Stub, t::Float64)
    n = floor(Int, (Float64(t) - s.clock) / s.dt)
    n <= 0 && return nothing
    s.clock += n * s.dt
    s.nsteps += n
    # a real step takes time; make any concurrent read observable
    sleep(0.002)
    return nothing
end

#####
##### ITime variants
#####
# Every real config sets `use_itime: true`, because Float comparisons of time
# drift over long runs (observed with Float32 time). Under ITime the
# Oceananigans clock holds a `DateTime`, so any code that coerces the coupler
# time to Float64 before calling `will_step` dispatches to the wrong method and
# fails with `Float64 - DateTime`. These stubs deliberately provide ONLY ITime
# methods, so such a coercion throws here rather than in a 45-minute GPU run.

mutable struct StubOceanIT <: Interfacer.AbstractOceanSimulation
    clock::Dates.DateTime
    dt::Dates.Second
    epoch::Dates.DateTime
    nsteps::Int
end
mutable struct StubIceIT <: Interfacer.AbstractSeaIceSimulation
    clock::Dates.DateTime
    dt::Dates.Second
    epoch::Dates.DateTime
    nsteps::Int
end
const StubIT = Union{StubOceanIT, StubIceIT}

Interfacer.sim_dt(s::StubIT) = Float64(Dates.value(s.dt))
Interfacer.will_step(s::StubIT, t::ITime) = (date(t) - s.clock) >= s.dt
function Interfacer.step!(s::StubIT, t::ITime)
    n = Dates.value(date(t) - s.clock) ÷ (1000 * Dates.value(s.dt))
    n <= 0 && return nothing
    s.clock += n * s.dt
    s.nsteps += n
    sleep(0.002)
    return nothing
end

function build_cs_itime(; dt_cpl, dt_slow, overlap, prime = false)
    epoch = Dates.DateTime(2010, 1, 1)
    ocean = StubOceanIT(epoch, Dates.Second(Int(dt_slow)), epoch, 0)
    ice = StubIceIT(epoch, Dates.Second(Int(dt_slow)), epoch, 0)
    Δt = ITime(Int64(dt_cpl), period = Dates.Second(1), epoch = epoch)
    t0 = ITime(Int64(0), period = Dates.Second(1), epoch = epoch)
    cs = Interfacer.CoupledSimulation{Float64}(
        epoch, nothing, nothing, (t0, t0), Δt, Ref(t0), Ref(0), Ref(-1),
        (; ice_sim = ice, ocean_sim = ocean),
        (), (;), nothing, nothing, false, true, overlap, prime,
        Ref{Any}(nothing), Ref{Any}(nothing), (;),
    )
    return cs, ocean, ice
end

function build_cs(; dt_cpl, dt_slow, overlap, prime = false)
    ocean = StubOcean(0.0, dt_slow, 0)
    ice = StubIce(0.0, dt_slow, 0)
    cs = Interfacer.CoupledSimulation{Float64}(
        nothing,                       # start_date
        nothing,                       # fields
        nothing,                       # conservation_checks
        (0.0, Inf),                    # tspan
        dt_cpl,                        # Δt_cpl
        Ref(0.0),                      # t
        Ref(0),                        # step
        Ref(-1),                       # prev_checkpoint_t
        (; ice_sim = ice, ocean_sim = ocean),
        (),                            # callbacks
        (;),                           # dir_paths
        nothing,                       # thermo_params
        nothing,                       # diags_handler
        false,                         # save_cache
        true,                          # step_concurrently
        overlap,                       # overlap_slow_surfaces
        prime,                         # prime_slow_surfaces
        Ref{Any}(nothing),             # slow_task
        Ref{Any}(nothing),             # slow_progress
        (;),                           # flux_accumulators
    )
    return cs, ocean, ice
end

"Mirror the scheduling in SimCoordinator.step! (the parts that touch slow sims)."
function drive!(cs, ocean, nsteps)
    log = NamedTuple[]
    for _ in 1:nsteps
        cs.t[] += cs.Δt_cpl
        cs.step[] += 1

        due() = cs.prime_slow_surfaces ? SimCoordinator.at_slow_boundary(cs) :
                SimCoordinator.slow_surfaces_due(cs)
        if cs.overlap_slow_surfaces && due()
            FieldExchanger.wait_slow_sims!(cs)
        end

        frozen = FieldExchanger.slow_step_in_flight(cs)

        # step_model_sims! advances the slow group only when neither frozen nor overlapping
        clock_before = ocean.clock
        (frozen || cs.overlap_slow_surfaces) ||
            FieldExchanger.step_slow_sims!(cs.model_sims, cs.t[])
        stepped_sync = ocean.clock != clock_before
        # Sample where `exchange!` runs, i.e. after the stepping phase. In the
        # baseline the slow sims step inside the loop, so sampling before that
        # would understate what the atmosphere actually sees.
        clock_at_exchange = ocean.clock

        launched = false
        if cs.overlap_slow_surfaces && !frozen && due()
            FieldExchanger.launch_slow_sims!(cs)
            launched = true
        end
        push!(
            log,
            (; step = cs.step[], t = cs.t[], frozen, stepped_sync, launched,
               ocean_clock_seen = clock_at_exchange),
        )
    end
    FieldExchanger.wait_slow_sims!(cs)
    return log
end

function report(label; dt_cpl, dt_slow, nsteps, overlap)
    cs, ocean, ice = build_cs(; dt_cpl, dt_slow, overlap)
    local log
    try
        log = drive!(cs, ocean, nsteps)
    catch e
        println("\n### $label -> THREW: ", e)
        return false
    end
    k = Int(dt_slow / dt_cpl)
    println("\n### $label   (dt_cpl=$dt_cpl, dt_slow=$dt_slow, k=$k, overlap=$overlap)")
    println("  step |    t | frozen | sync-stepped | launched | ocean clock at entry")
    for r in log
        println("  ", lpad(r.step, 4), " | ", lpad(Int(r.t), 4), " | ",
                lpad(r.frozen, 6), " | ", lpad(r.stepped_sync, 12), " | ",
                lpad(r.launched, 8), " | ", lpad(Int(r.ocean_clock_seen), 6))
    end

    # invariants
    ok = true
    expected_steps = Int(floor(nsteps * dt_cpl / dt_slow))
    if ocean.nsteps != expected_steps
        println("  FAIL: ocean took $(ocean.nsteps) steps, expected $expected_steps")
        ok = false
    end
    if ice.nsteps != ocean.nsteps
        println("  FAIL: ice took $(ice.nsteps) steps, ocean took $(ocean.nsteps)")
        ok = false
    end
    if any(r -> r.frozen && r.stepped_sync, log)
        println("  FAIL: a slow sim was stepped synchronously while a step was in flight")
        ok = false
    end
    if ocean.clock > nsteps * dt_cpl
        println("  FAIL: ocean ran past the coupler ($(ocean.clock) > $(nsteps*dt_cpl))")
        ok = false
    end
    overlapped = count(r -> r.frozen, log)
    println("  ocean steps: $(ocean.nsteps) (expected $expected_steps), ",
            "final clock: $(ocean.clock), coupler: $(nsteps*dt_cpl)")
    println("  coupling steps overlapped with an in-flight slow step: $overlapped")
    println(ok ? "  PASS" : "  FAIL")
    return ok
end

"Same schedule, but with ITime times and stubs that reject Float64 coercion."
function report_itime(label; dt_cpl, dt_slow, nsteps, overlap)
    cs, ocean, ice = build_cs_itime(; dt_cpl, dt_slow, overlap)
    k = Int(dt_slow / dt_cpl)
    try
        drive!(cs, ocean, nsteps)
    catch e
        println("\n### $label (ITime)  -> THREW: ", sprint(showerror, e))
        return false
    end
    expected = Int(floor(nsteps * dt_cpl / dt_slow))
    ok = ocean.nsteps == expected && ice.nsteps == ocean.nsteps
    println("\n### $label (ITime)   (dt_cpl=$dt_cpl, dt_slow=$dt_slow, k=$k, overlap=$overlap)")
    println("  ocean steps: $(ocean.nsteps) (expected $expected), ice: $(ice.nsteps)")
    println(ok ? "  PASS" : "  FAIL")
    return ok
end

"""
Compare the ocean state the atmosphere would see, step by step, across the three
schedules. Priming claims to remove the extra lag that plain overlap introduces;
this checks that claim directly rather than inferring it.
"""
function report_lag(; dt_cpl, dt_slow, nsteps)
    k = Int(dt_slow / dt_cpl)
    seen = Dict{String, Vector{Float64}}()
    counts = Dict{String, Int}()
    for (label, overlap, prime) in
        (("baseline", false, false), ("overlap", true, false), ("primed", true, true))
        cs, ocean, _ = build_cs(; dt_cpl, dt_slow, overlap, prime)
        if prime
            # what the constructor does when prime_slow_surfaces is set
            FieldExchanger.step_slow_sims!(cs.model_sims, cs.t[] + k * cs.Δt_cpl)
        end
        log = drive!(cs, ocean, nsteps)
        seen[label] = [r.ocean_clock_seen for r in log]
        counts[label] = ocean.nsteps
    end

    println("\n### ocean time the atmosphere sees, by coupling step (k=$k)")
    println("  step |  t_n | baseline | overlap | primed")
    for i in 1:nsteps
        println("  ", lpad(i, 4), " | ", lpad(Int(i * dt_cpl), 4), " | ",
                lpad(Int(seen["baseline"][i]), 8), " | ", lpad(Int(seen["overlap"][i]), 7),
                " | ", lpad(Int(seen["primed"][i]), 6))
    end
    # after the first window, primed should match baseline exactly
    tail = (k + 1):nsteps
    primed_ok = all(seen["primed"][i] == seen["baseline"][i] for i in tail)
    overlap_lags = any(seen["overlap"][i] < seen["baseline"][i] for i in tail)
    println("  ocean steps taken: ", [l => counts[l] for l in ("baseline", "overlap", "primed")])
    println("  primed matches baseline after the first window: ", primed_ok)
    println("  plain overlap lags baseline: ", overlap_lags)
    ok = primed_ok && overlap_lags
    println(ok ? "  PASS" : "  FAIL")
    return ok
end

allok = true
allok &= report("k=1  (dt_ocean == dt_cpl)"; dt_cpl = 360.0, dt_slow = 360.0, nsteps = 8, overlap = true)
allok &= report("k=5  (the usual case)";     dt_cpl = 360.0, dt_slow = 1800.0, nsteps = 20, overlap = true)
allok &= report("k=2";                        dt_cpl = 360.0, dt_slow = 720.0, nsteps = 10, overlap = true)
allok &= report("k=1, overlap OFF (control)"; dt_cpl = 360.0, dt_slow = 360.0, nsteps = 8, overlap = false)
allok &= report_itime("k=1  (dt_ocean == dt_cpl)"; dt_cpl = 360.0, dt_slow = 360.0, nsteps = 8, overlap = true)
allok &= report_itime("k=5  (the usual case)";     dt_cpl = 360.0, dt_slow = 1800.0, nsteps = 20, overlap = true)
allok &= report_lag(; dt_cpl = 360.0, dt_slow = 1800.0, nsteps = 15)
println("\n", allok ? "ALL PASS" : "FAILURES ABOVE")
