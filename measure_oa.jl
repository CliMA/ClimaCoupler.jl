# Measure O/A: the cost of one ocean/sea-ice step against one coupling step of
# atmos+land work. That ratio decides what overlapping can buy, since overlap
# hides at most one window of atmos work (k*A) regardless of how large O is:
#
#   saving = min(k*A, O) / (k*A + O)
#
# which peaks when O ~= k*A and falls away in both directions.
#
#   julia --project=experiments/CMIP -t 4 measure_oa.jl [config_name]
#
# Method notes -- the previous version of this script got all three wrong:
#
#  1. Timing. CUDA.@elapsed (what ClimaComms.@elapsed uses on a CUDADevice)
#     brackets GPU stream work with events. It can miss CPU-side time between
#     launches, and it cannot see work pushed onto other streams by the
#     concurrent path. Here we take wall-clock time with an explicit device
#     synchronize on both sides, which captures everything.
#
#  2. Warmup. With k = 5 the ocean does not step during the first four coupling
#     steps, so a 4-step warmup left its kernels uncompiled and the first ocean
#     sample was 388 s of JIT. Warm up past 2k steps so the ocean has stepped
#     at least twice.
#
#  3. Decomposition. Timing component groups out of band misses coupler
#     overhead and cannot be reconciled with the run's real per-step cost.
#     Instead we time whole coupling steps through the real code path and take
#     O as the difference between steps where the ocean does and does not step.
#     That requires the SEQUENTIAL path (step_concurrently: false), or the two
#     groups overlap and the difference measures max(A,O) rather than O.
#
# The script cross-checks its own numbers against the observed per-coupling-step
# wall time and says so if they do not reconcile.

include(joinpath(@__DIR__, "experiments", "CMIP", "code_loading.jl"))
import Dates
import CUDA
import Statistics: median

const Interfacer = ClimaCoupler.Interfacer
const FieldExchanger = ClimaCoupler.FieldExchanger
const FluxCalculator = ClimaCoupler.FluxCalculator
const SimCoordinator = ClimaCoupler.SimCoordinator

const NAME = isempty(ARGS) ? "oa_serial" : ARGS[1]
const CONFIG = joinpath("config", "ci_configs", "$(NAME).yml")

empty!(ARGS)
append!(ARGS, ["--config_file", CONFIG])

"Wall-clock seconds with the device drained before and after."
function timed(f)
    CUDA.synchronize()
    t0 = time_ns()
    f()
    CUDA.synchronize()
    return (time_ns() - t0) / 1e9
end

@info "building $NAME"
cs = Interfacer.CoupledSimulation(CONFIG)
sims = cs.model_sims

Δt_cpl = Float64(float(cs.Δt_cpl))
dt_slow = FieldExchanger.slow_sim_dt(cs)
k = round(Int, dt_slow / Δt_cpl)
@info "timesteps" Δt_cpl dt_slow k step_concurrently = cs.step_concurrently

if cs.step_concurrently
    @warn "step_concurrently is TRUE: the ocean overlaps the atmos group, so the " *
          "per-step difference measures max(A,O), not O. Use a sequential config."
end

"True if the slow surfaces will step when the coupler advances to the next time."
slow_steps_next(cs) = any(
    sim -> Interfacer.is_overlapped(sim) && Interfacer.will_step(sim, cs.t[] + cs.Δt_cpl),
    values(cs.model_sims),
)

# --- warmup: past 2k coupling steps so the ocean has stepped twice and every
# --- kernel on both paths is compiled
n_warm = 2k + 2
@info "warming up ($n_warm coupling steps, so the ocean steps twice)"
for i in 1:n_warm
    dt = timed(() -> SimCoordinator.step!(cs))
    @info "warmup $i/$n_warm: $(round(dt, digits = 3)) s"
end

# --- measure whole coupling steps through the real path ---
n_meas = 3k
@info "measuring ($n_meas coupling steps)"
with_ocean = Float64[]
without_ocean = Float64[]
for i in 1:n_meas
    will = slow_steps_next(cs)
    dt = timed(() -> SimCoordinator.step!(cs))
    push!(will ? with_ocean : without_ocean, dt)
    @info "step $i: $(round(dt, digits = 3)) s   ocean stepped: $will"
end

# --- direct group timings, as a cross-check on the difference ---
t_atmos_group = timed(
    () -> Interfacer.step!(
        sims.land_sim,
        sims.atmos_sim,
        cs.t[],
        cs.fields,
        cs.thermo_params,
    ),
)
t_exchange = timed(function ()
    FieldExchanger.update_surface_fractions!(cs)
    FieldExchanger.exchange!(cs)
    FluxCalculator.turbulent_fluxes!(cs)
    FluxCalculator.ocean_seaice_fluxes!(cs)
end)

step_no_ocean = isempty(without_ocean) ? NaN : median(without_ocean)
step_ocean = isempty(with_ocean) ? NaN : median(with_ocean)
O = step_ocean - step_no_ocean
A = step_no_ocean          # a coupling step with no ocean work: atmos+land+coupler
kA = k * A

println("\n", "="^70)
println("config: $NAME    k = $k    step_concurrently = $(cs.step_concurrently)")
println("  coupling step, ocean idle    : ", round(step_no_ocean, digits = 3), " s   n=",
        length(without_ocean))
println("  coupling step, ocean steps   : ", round(step_ocean, digits = 3), " s   n=",
        length(with_ocean))
println("  => O (one ocean+ice step)    : ", round(O, digits = 3), " s")
println("  => A (coupling step w/o ocean): ", round(A, digits = 3), " s")
println("-"^70)
println("  cross-check, timed directly:")
println("    atmos+land group           : ", round(t_atmos_group, digits = 3), " s")
println("    exchange + fluxes          : ", round(t_exchange, digits = 3), " s")
println("    sum vs 'ocean idle' step   : ",
        round(t_atmos_group + t_exchange, digits = 3), " s vs ",
        round(step_no_ocean, digits = 3), " s")
println("-"^70)
println("  O/A = ", round(O / A, digits = 3))
println("  k*A = ", round(kA, digits = 3), " s")
saving = O > 0 ? min(kA, O) / (kA + O) : 0.0
println("  predicted overlap saving = min(kA,O)/(kA+O) = ", round(100 * saving, digits = 1), "%")
println("  (peaks at 50% when O == k*A)")
println("-"^70)
mean_step = (sum(with_ocean) + sum(without_ocean)) / n_meas
println("  mean coupling step here      : ", round(mean_step, digits = 3), " s")
println("  reconcile: A + O/k           = ", round(A + O / k, digits = 3), " s")
if !isfinite(O) || abs((A + O / k) - mean_step) > 0.15 * mean_step
    println("  *** THESE DO NOT RECONCILE -- do not trust O/A above ***")
else
    println("  reconciles within 15% -- decomposition is consistent")
end
println("="^70)
