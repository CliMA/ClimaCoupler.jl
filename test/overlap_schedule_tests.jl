# Exercise the overlapped slow-surface schedule without running a model.
#
# The launch/join logic assumes a window: launch at the end of one coupling step,
# join at the top of a later one. When dt_ocean == dt_cpl (k = 1) those collapse
# onto consecutive steps, which is legal and needs checking.
#
#   julia --project=experiments/CMIP test_overlap_schedule.jl

using Test
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

# Overlapping is opted into per concrete type, so stubs must say so too --
# prescribed and slab surfaces deliberately do not.
Interfacer.is_overlapped(::Stub) = true

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

Interfacer.is_overlapped(::StubIT) = true

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
        epoch,
        nothing,
        nothing,
        (t0, t0),
        Δt,
        Ref(t0),
        Ref(0),
        Ref(-1),
        (; ice_sim = ice, ocean_sim = ocean),
        (),
        (;),
        nothing,
        nothing,
        false,
        true,
        overlap,
        false,               # prime_fast_group
        Ref{Any}(nothing),   # slow_task
        Ref{Any}(nothing),   # slow_progress
        (;),
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
        false,                       # prime_fast_group
        Ref{Any}(nothing),             # slow_task
        Ref{Any}(nothing),             # slow_progress
        (;),                           # flux_accumulators
    )
    return cs, ocean, ice
end

"""
Drive the schedule using the very functions `SimCoordinator.step!` uses, rather
than restating their order here. An earlier version of this file duplicated that
order and went stale the moment `step!` changed, failing for a reason unrelated
to the code under test.

The rest of `step!` (exchange, fluxes, diagnostics) needs a real
`CoupledSimulation`, so it is not driven here; only the slow-surface schedule is.
"""
function drive!(cs, ocean, nsteps; prime_steps::Int = 0)
    log = NamedTuple[]
    # The ocean state the coupler is *permitted* to read. While a step is in
    # flight the coupler is frozen out, so reading ocean.clock directly would
    # race the task and report a value the real coupler could not have used.
    visible = ocean.clock
    for i in 1:nsteps
        cs.t[] += cs.Δt_cpl
        cs.step[] += 1

        frozen = SimCoordinator.join_slow_if_due!(cs)

        clock_before = ocean.clock
        SimCoordinator.skip_slow_stepping(cs, frozen) ||
            FieldExchanger.step_slow_sims!(cs.model_sims, cs.t[])
        stepped_sync = ocean.clock != clock_before

        # Safe to read whenever nothing is in flight; covers both the baseline
        # (stepped in the loop just above) and a boundary step (just joined).
        frozen || (visible = ocean.clock)

        # Mirrors `step!(cs; suppress_slow_launch = true)`: during priming the
        # fast group advances but the slow group is not started.
        launched = i <= prime_steps ? false : SimCoordinator.launch_slow_if_due!(cs, frozen)
        push!(
            log,
            (;
                step = cs.step[],
                t = cs.t[],
                frozen,
                stepped_sync,
                launched,
                ocean_clock_seen = visible,
            ),
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
        println(
            "  ",
            lpad(r.step, 4),
            " | ",
            lpad(Int(r.t), 4),
            " | ",
            lpad(r.frozen, 6),
            " | ",
            lpad(r.stepped_sync, 12),
            " | ",
            lpad(r.launched, 8),
            " | ",
            lpad(Int(r.ocean_clock_seen), 6),
        )
    end

    expected_steps = Int(floor(nsteps * dt_cpl / dt_slow))
    @test ocean.nsteps == expected_steps
    @test ice.nsteps == ocean.nsteps
    # nothing may be stepped in the loop while its task owns it
    @test !any(r -> r.frozen && r.stepped_sync, log)
    # the slow group must never run past the coupler
    @test ocean.clock <= nsteps * dt_cpl
    overlapped = count(r -> r.frozen, log)
    println(
        "  ocean steps: $(ocean.nsteps) (expected $expected_steps), ",
        "final clock: $(ocean.clock), coupler: $(nsteps*dt_cpl)",
    )
    println("  coupling steps overlapped with an in-flight slow step: $overlapped")
    return nothing
end

"Same schedule, but with ITime times and stubs that reject Float64 coercion."
function report_itime(label; dt_cpl, dt_slow, nsteps, overlap)
    cs, ocean, ice = build_cs_itime(; dt_cpl, dt_slow, overlap)
    k = Int(dt_slow / dt_cpl)
    # Any coercion of coupler time to Float64 throws here, because the ITime
    # stubs deliberately define only ITime methods.
    @test (drive!(cs, ocean, nsteps); true)
    expected = Int(floor(nsteps * dt_cpl / dt_slow))
    println(
        "\n### $label (ITime)   (dt_cpl=$dt_cpl, dt_slow=$dt_slow, k=$k, overlap=$overlap)",
    )
    println("  ocean steps: $(ocean.nsteps) (expected $expected), ice: $(ice.nsteps)")
    @test ocean.nsteps == expected
    @test ice.nsteps == ocean.nsteps
    return nothing
end

"""
Compare the ocean state the atmosphere would see, step by step, with and without
overlapping, so the extra lag overlapping introduces is shown rather than
inferred.
"""
function report_lag(; dt_cpl, dt_slow, nsteps)
    k = Int(dt_slow / dt_cpl)
    seen = Dict{String, Vector{Float64}}()
    counts = Dict{String, Int}()
    for (label, overlap) in (("baseline", false), ("overlap", true))
        cs, ocean, _ = build_cs(; dt_cpl, dt_slow, overlap)
        log = drive!(cs, ocean, nsteps)
        seen[label] = [r.ocean_clock_seen for r in log]
        counts[label] = ocean.nsteps
    end

    println("\n### ocean time the atmosphere sees, by coupling step (k=$k)")
    println("  step |  t_n | baseline | overlap")
    for i in 1:nsteps
        println(
            "  ",
            lpad(i, 4),
            " | ",
            lpad(Int(i * dt_cpl), 4),
            " | ",
            lpad(Int(seen["baseline"][i]), 8),
            " | ",
            lpad(Int(seen["overlap"][i]), 7),
        )
    end
    tail = (k + 1):nsteps
    overlap_lags = any(seen["overlap"][i] < seen["baseline"][i] for i in tail)
    println("  ocean steps taken: ", [l => counts[l] for l in ("baseline", "overlap")])
    println("  plain overlap lags baseline: ", overlap_lags)
    @test overlap_lags
    return nothing
end

report(
    "k=1  (dt_ocean == dt_cpl)";
    dt_cpl = 360.0,
    dt_slow = 360.0,
    nsteps = 8,
    overlap = true,
)
report(
    "k=5  (the usual case)";
    dt_cpl = 360.0,
    dt_slow = 1800.0,
    nsteps = 20,
    overlap = true,
)
report("k=2"; dt_cpl = 360.0, dt_slow = 720.0, nsteps = 10, overlap = true)
report(
    "k=1, overlap OFF (control)";
    dt_cpl = 360.0,
    dt_slow = 360.0,
    nsteps = 8,
    overlap = false,
)
report_itime(
    "k=1  (dt_ocean == dt_cpl)";
    dt_cpl = 360.0,
    dt_slow = 360.0,
    nsteps = 8,
    overlap = true,
)
report_itime(
    "k=5  (the usual case)";
    dt_cpl = 360.0,
    dt_slow = 1800.0,
    nsteps = 20,
    overlap = true,
)
"""
Restart scenario. `cs.step[]` always begins again at zero on a restart, and
`checkpoint_sims` joins any in-flight slow step before saving, so a checkpoint
taken mid-window leaves the slow group an arbitrary amount ahead of the coupler
-- not a whole window. The schedule therefore cannot be a count of coupling
steps; it has to be re-derived from where the components actually are.
"""
report_lag(; dt_cpl = 360.0, dt_slow = 1800.0, nsteps = 15)
# clean boundary checkpoint: slow group exactly one window ahead
# mid-window checkpoint: slow group ahead by an amount that is NOT a whole window

# ---------------------------------------------------------------------------
# While a slow step is in flight, the coupler must not compute turbulent fluxes
# for the overlapped surfaces: doing so reads surface state the task is writing,
# and for a surface with no accumulator it writes fluxes back into the model.
# Their share is parked at the last sync and restored instead.
#
# `csf` here is a NamedTuple of plain vectors: `turbulent_fluxes!` only needs
# `propertynames`, `getproperty`, `fill!` and broadcast, so a full ClimaCore
# space is unnecessary and would obscure what is being checked.

const FluxCalculator = ClimaCoupler.FluxCalculator

mutable struct CountingAtmos <: Interfacer.AbstractAtmosSimulation end

mutable struct CountingSurface <: Interfacer.AbstractSurfaceSimulation
    overlapped::Bool
    contribution::Float64
    ncalls::Int
end
CountingSurface(o, c) = CountingSurface(o, c, 0)
Interfacer.is_overlapped(s::CountingSurface) = s.overlapped

function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::CountingSurface,
    ::Interfacer.AbstractAtmosSimulation,
    thermo_params,
    accumulator = nothing,
)
    sim.ncalls += 1
    csf.F_sh .+= sim.contribution
    return nothing
end

@testset "overlapped surfaces are skipped while frozen" begin
    names = (
        :F_turb_ρτxz,
        :F_turb_ρτyz,
        :F_lh,
        :F_sh,
        :F_turb_moisture,
        :slow_F_turb_ρτxz,
        :slow_F_turb_ρτyz,
        :slow_F_lh,
        :slow_F_sh,
        :slow_F_turb_moisture,
    )
    csf = NamedTuple{names}(Tuple(zeros(3) for _ in names))

    ocean = CountingSurface(true, 2.0)     # overlapped
    land = CountingSurface(false, 5.0)     # fast
    sims = (; atmos_sim = CountingAtmos(), ocean_sim = ocean, land_sim = land)

    # Sync step: both surfaces contribute, and the overlapped share is parked.
    FluxCalculator.turbulent_fluxes!(csf, sims, nothing; slow_frozen = false)
    @test ocean.ncalls == 1
    @test land.ncalls == 1
    @test all(csf.F_sh .== 7.0)
    @test all(csf.slow_F_sh .== 2.0)

    # Frozen step: the overlapped surface must not be touched, and its parked
    # share must still reach the total.
    FluxCalculator.turbulent_fluxes!(csf, sims, nothing; slow_frozen = true)
    @test ocean.ncalls == 1          # not called again
    @test land.ncalls == 2           # fast surfaces still computed live
    @test all(csf.F_sh .== 7.0)      # parked slow share restored, not dropped
    @test all(csf.slow_F_sh .== 2.0) # cache untouched while frozen

    # Back in sync: the overlapped surface is read again and the cache refreshed.
    ocean.contribution = 3.0
    FluxCalculator.turbulent_fluxes!(csf, sims, nothing; slow_frozen = false)
    @test ocean.ncalls == 2
    @test all(csf.F_sh .== 8.0)
    @test all(csf.slow_F_sh .== 3.0)
end

@testset "is_overlapped is opt-in by concrete type" begin
    # Prescribed and slab surfaces must not be dragged into the overlapped group:
    # they are cheap, so they would pay a coupling step of lag for nothing.
    @test !Interfacer.is_overlapped(CountingSurface(false, 0.0))
    @test Interfacer.is_overlapped(CountingSurface(true, 0.0))
    @test !Interfacer.is_overlapped(
        ClimaCoupler.Interfacer.SurfaceStub((; area_fraction = nothing)),
    )
end


# ---------------------------------------------------------------------------
# `prime_fast_group` runs the fast group k-1 coupling steps before the slow
# group starts, so a slow step is launched carrying a whole window of forcing.
# The cost is that the two groups' model times stay offset by k-1 steps.

@testset "priming the fast group delays the first slow launch" begin
    dt_cpl, dt_slow = 360.0, 1800.0
    k = Int(dt_slow / dt_cpl)

    cs_p, ocean_p, _ = build_cs(; dt_cpl, dt_slow, overlap = true)
    log_p = drive!(cs_p, ocean_p, 2k; prime_steps = k - 1)

    cs_n, ocean_n, _ = build_cs(; dt_cpl, dt_slow, overlap = true)
    log_n = drive!(cs_n, ocean_n, 2k)

    first_launch(log) = findfirst(r -> r.launched, log)
    println("\n### fast-group priming (k=$k)")
    println("  first slow launch, unprimed: step ", first_launch(log_n))
    println("  first slow launch, primed:   step ", first_launch(log_p))
    println(
        "  ocean steps over $(2k) coupling steps: ",
        "unprimed ",
        ocean_n.nsteps,
        ", primed ",
        ocean_p.nsteps,
    )

    # Nothing is launched while the fast group is being primed.
    @test all(!r.launched for r in log_p[1:(k - 1)])
    # And the slow group still runs: it is delayed, not suppressed.
    @test first_launch(log_p) !== nothing
    @test first_launch(log_p) >= k
    # Priming costs the slow group at most one step over the window shown.
    @test ocean_p.nsteps >= ocean_n.nsteps - 1
end
