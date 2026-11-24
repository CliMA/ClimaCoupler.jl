# The overlapped slow-surface schedule, driven with stub ocean and sea ice models.

using Test
using ClimaCoupler
import ClimaComms
import ClimaCore as CC
import Dates
import ClimaUtilities.TimeManager: ITime, date
const Interfacer = ClimaCoupler.Interfacer
const FieldExchanger = ClimaCoupler.FieldExchanger
const FluxCalculator = ClimaCoupler.FluxCalculator
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

Interfacer.is_overlapped(::Stub) = true
Interfacer.sim_dt(s::Stub) = s.dt
Interfacer.will_step(s::Stub, t::Float64) = (t - s.clock) >= s.dt
function Interfacer.step!(s::Stub, t::Float64)
    n = floor(Int, (t - s.clock) / s.dt)
    n <= 0 && return nothing
    s.clock += n * s.dt
    s.nsteps += n
    # A step that takes time makes a concurrent read observable.
    sleep(0.002)
    return nothing
end

# The ITime stubs define only ITime methods, so any coercion of coupler time to
# Float64 before `will_step` or `step!` throws here.
mutable struct StubOceanIT <: Interfacer.AbstractOceanSimulation
    clock::Dates.DateTime
    dt::Dates.Second
    nsteps::Int
end
mutable struct StubIceIT <: Interfacer.AbstractSeaIceSimulation
    clock::Dates.DateTime
    dt::Dates.Second
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

# Coupler fields, read only for the device the slow task synchronizes.
stub_fields() = CC.Fields.zeros(
    CC.CommonSpaces.CubedSphereSpace(
        Float64;
        radius = 6.371e6,
        n_quad_points = 2,
        h_elem = 1,
    ),
)

function build_cs(ocean, ice, t0, Δt_cpl, t_end, overlap)
    return Interfacer.CoupledSimulation{Float64}(
        nothing,                       # start_date
        stub_fields(),                 # fields
        nothing,                       # conservation_checks
        (t0, t_end),                   # tspan
        Δt_cpl,                        # Δt_cpl
        Ref(t0),                       # t
        Ref(0),                        # step
        Ref(-1),                       # prev_checkpoint_t
        (; ice_sim = ice, ocean_sim = ocean),
        (),                            # callbacks
        (;),                           # dir_paths
        nothing,                       # thermo_params
        nothing,                       # diags_handler
        false,                         # save_cache
        Interfacer.SteppingSchedule(true, overlap),
        (;),                           # flux_accumulators
    )
end

function build_cs(; k, nsteps, overlap, dt_cpl = 360.0)
    ocean, ice = StubOcean(0.0, k * dt_cpl, 0), StubIce(0.0, k * dt_cpl, 0)
    return build_cs(ocean, ice, 0.0, dt_cpl, nsteps * dt_cpl, overlap), ocean, ice
end

function build_cs_itime(; k, nsteps, overlap, dt_cpl = 360)
    epoch = Dates.DateTime(2010, 1, 1)
    ocean = StubOceanIT(epoch, Dates.Second(k * dt_cpl), 0)
    ice = StubIceIT(epoch, Dates.Second(k * dt_cpl), 0)
    t(seconds) = ITime(seconds, period = Dates.Second(1), epoch = epoch)
    return build_cs(ocean, ice, t(0), t(dt_cpl), t(nsteps * dt_cpl), overlap), ocean, ice
end

"""
Run the slow-surface schedule for `nsteps` coupling steps, and return the ocean
clock the coupler may read after each one.
"""
function drive!(cs, ocean, nsteps)
    visible = ocean.clock
    seen = typeof(visible)[]
    for _ in 1:nsteps
        cs.t[] += cs.Δt_cpl
        cs.step[] += 1
        cs.schedule.overlap_slow ||
            FieldExchanger.step_slow_sims!(cs.model_sims, cs.t[], ClimaComms.device(cs))
        frozen = SimCoordinator.join_slow_if_due!(cs)
        frozen || (visible = ocean.clock)
        push!(seen, visible)
        SimCoordinator.launch_slow_if_due!(cs, frozen)
    end
    FieldExchanger.wait_slow_sims!(cs)
    return seen
end

@testset "one slow step per window, never past the end of the run" begin
    # nsteps = 19 ends the run part-way through a window.
    for (k, nsteps) in ((1, 8), (2, 10), (5, 20), (5, 19)), overlap in (false, true)
        cs, ocean, ice = build_cs(; k, nsteps, overlap)
        drive!(cs, ocean, nsteps)
        @test ocean.nsteps == nsteps ÷ k
        @test ice.nsteps == ocean.nsteps
        @test ocean.clock <= nsteps * cs.Δt_cpl
    end
end

@testset "one slow step per window with ITime" begin
    for k in (1, 5)
        cs, ocean, ice = build_cs_itime(; k, nsteps = 20, overlap = true)
        drive!(cs, ocean, 20)
        @test ocean.nsteps == 20 ÷ k
        @test ice.nsteps == ocean.nsteps
    end
end

@testset "overlap delays the ocean the coupler sees by max(k - 1, 1) steps" begin
    for k in (1, 2, 5)
        nsteps = 4k
        sequential, overlapped = map((false, true)) do overlap
            cs, ocean, _ = build_cs(; k, nsteps, overlap)
            drive!(cs, ocean, nsteps)
        end
        lag = max(k - 1, 1)
        @test overlapped[(lag + 1):end] == sequential[1:(end - lag)]
    end
end

# Under overlap the coupler must not compute turbulent fluxes for an overlapped
# surface while its step is in flight; the share cached at the last join is restored.
# `csf` is a NamedTuple of vectors, which is all `turbulent_fluxes!` needs.

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
    names = (:F_turb_ρτxz, :F_turb_ρτyz, :F_lh, :F_sh, :F_turb_moisture)
    names = (names..., map(Interfacer.slow_cache_name, names)...)
    csf = NamedTuple{names}(Tuple(zeros(3) for _ in names))

    ocean = CountingSurface(true, 2.0)
    land = CountingSurface(false, 5.0)
    sims = (; atmos_sim = CountingAtmos(), ocean_sim = ocean, land_sim = land)

    # Both surfaces contribute, and the overlapped share is cached.
    FluxCalculator.turbulent_fluxes!(csf, sims, nothing; slow_frozen = false)
    @test (ocean.ncalls, land.ncalls) == (1, 1)
    @test all(csf.F_sh .== 7.0)

    # Frozen: the overlapped surface is not called, and its cached share still
    # reaches the total.
    FluxCalculator.turbulent_fluxes!(csf, sims, nothing; slow_frozen = true)
    @test (ocean.ncalls, land.ncalls) == (1, 2)
    @test all(csf.F_sh .== 7.0)

    # After the join the overlapped surface is read again.
    ocean.contribution = 3.0
    FluxCalculator.turbulent_fluxes!(csf, sims, nothing; slow_frozen = false)
    @test ocean.ncalls == 2
    @test all(csf.F_sh .== 8.0)
end
