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

function build_cs_itime(; dt_cpl, dt_slow, overlap)
    epoch = Dates.DateTime(2010, 1, 1)
    ocean = StubOceanIT(epoch, Dates.Second(Int(dt_slow)), epoch, 0)
    ice = StubIceIT(epoch, Dates.Second(Int(dt_slow)), epoch, 0)
    Δt = ITime(Int64(dt_cpl), period = Dates.Second(1), epoch = epoch)
    t0 = ITime(Int64(0), period = Dates.Second(1), epoch = epoch)
    cs = Interfacer.CoupledSimulation{Float64}(
        epoch, nothing, nothing, (t0, t0), Δt, Ref(t0), Ref(0), Ref(-1),
        (; ice_sim = ice, ocean_sim = ocean),
        (), (;), nothing, nothing, false, true, overlap,
        Ref{Any}(nothing), Ref{Any}(nothing), (;),
    )
    return cs, ocean, ice
end

function build_cs(; dt_cpl, dt_slow, overlap)
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

        if cs.overlap_slow_surfaces && SimCoordinator.slow_surfaces_due(cs)
            FieldExchanger.wait_slow_sims!(cs)
        end
        frozen = FieldExchanger.slow_step_in_flight(cs)

        # step_model_sims! advances the slow group only when neither frozen nor overlapping
        clock_before = ocean.clock
        (frozen || cs.overlap_slow_surfaces) ||
            FieldExchanger.step_slow_sims!(cs.model_sims, cs.t[])
        stepped_sync = ocean.clock != clock_before

        launched = false
        if cs.overlap_slow_surfaces && !frozen && SimCoordinator.slow_surfaces_due(cs)
            FieldExchanger.launch_slow_sims!(cs)
            launched = true
        end
        push!(
            log,
            (; step = cs.step[], t = cs.t[], frozen, stepped_sync, launched,
               ocean_clock_seen = clock_before),
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

allok = true
allok &= report("k=1  (dt_ocean == dt_cpl)"; dt_cpl = 360.0, dt_slow = 360.0, nsteps = 8, overlap = true)
allok &= report("k=5  (the usual case)";     dt_cpl = 360.0, dt_slow = 1800.0, nsteps = 20, overlap = true)
allok &= report("k=2";                        dt_cpl = 360.0, dt_slow = 720.0, nsteps = 10, overlap = true)
allok &= report("k=1, overlap OFF (control)"; dt_cpl = 360.0, dt_slow = 360.0, nsteps = 8, overlap = false)
allok &= report_itime("k=1  (dt_ocean == dt_cpl)"; dt_cpl = 360.0, dt_slow = 360.0, nsteps = 8, overlap = true)
allok &= report_itime("k=5  (the usual case)";     dt_cpl = 360.0, dt_slow = 1800.0, nsteps = 20, overlap = true)
println("\n", allok ? "ALL PASS" : "FAILURES ABOVE")
