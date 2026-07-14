# experiments/CMIP/stability_diagnostics.jl
#
# Diagnoses the source of the ~3.6-day instability in the CMIP
# oceananigans+climaseaice configuration.
#
# Three tests run in sequence:
#
#   Test 1  – Run N_MONITOR coupling steps with a single CoupledSimulation,
#              printing field extrema each step.  Identifies which component
#              (ocean T/S/u, SIC, T_sfc, turbulent fluxes) first goes NaN or
#              unphysical.  Stops early at the first NaN.
#
#   Test 2A – Construct two CoupledSimulation objects back-to-back (before any
#              stepping) and compare their t=0 states.  All differences should
#              be exactly 0.0 if construction is deterministic and free of
#              global-state side-effects.
#
#   Test 2B – Construct a fresh CoupledSimulation *after* Test 1's run/crash
#              and compare its t=0 state against Test 1's t=0 state.
#              Nonzero differences → global-state pollution from the first run
#              (the two cs objects have different effective ICs, which explains
#              different trajectories).  All-zero differences → same ICs, so
#              the divergence must come from non-determinism inside step!.
#
# Configuration via environment variables (all optional):
#   STABILITY_CONFIG   path to YAML config file (default below)
#   STABILITY_N_STEPS  coupling steps for Test 1  (default 300 ≈ 3.6 days)
#   STABILITY_JOB_ID   job id used to name the output directory

include("code_loading.jl")

const CONFIG_FILE  = get(ENV, "STABILITY_CONFIG",
                         "config/ci_configs/cmip_oceananigans_climaseaice.yml")
const N_MONITOR    = parse(Int, get(ENV, "STABILITY_N_STEPS", "300"))
const JOB_ID       = get(ENV, "STABILITY_JOB_ID", "stability_diagnostics")
const ARTIFACT_DIR = joinpath("output", JOB_ID, "artifacts")
mkpath(ARTIFACT_DIR)

@info "Config      : $CONFIG_FILE"
@info "N_MONITOR   : $N_MONITOR"
@info "Artifact dir: $ARTIFACT_DIR"

# ─── helpers ──────────────────────────────────────────────────────────────────
#
# Oceananigans.interior / ClimaCore Fields may live on GPU; always materialize
# with Array(...) before extrema / NaN checks so diagnostics work on CPU and GPU.
# @sprintf requires a *literal* format string (no "a" * "b"); keep formats here.

# Materialize an Oceananigans interior slice to a CPU Vector.
function interior_vec(field, i, j, k)
    return vec(Array(Oceananigans.interior(field, i, j, k)))
end

# Materialize a ClimaCore coupler field to a CPU Vector.
function field_vec(field)
    return vec(Array(parent(field)))
end

function finite_extrema(x)
    xmin, xmax = extrema(x)
    return (Float64(xmin), Float64(xmax))
end

function has_nan(x)
    return any(!isfinite, x)
end

function snapshot_state(cs)
    ocean = cs.model_sims.ocean_sim.ocean.model
    ice   = cs.model_sims.ice_sim.ice.model
    Nz    = size(ocean.grid, 3)
    return (;
        T_oc      = interior_vec(ocean.tracers.T,       :, :, Nz),
        S_oc      = interior_vec(ocean.tracers.S,       :, :, Nz),
        u_oc      = interior_vec(ocean.velocities.u,    :, :, Nz),
        SIC       = interior_vec(ice.ice_concentration, :, :, 1),
        h_ice     = interior_vec(ice.ice_thickness,     :, :, 1),
        T_ice_top = interior_vec(ice.ice_thermodynamics.top_surface_temperature, :, :, 1),
        T_sfc     = field_vec(cs.fields.T_sfc),
        t         = Float64(cs.t[]),
    )
end

function write_comparison(io, label, s1, s2)
    println(io, "\n=== $label ===")
    all_zero = true
    for f in (:T_oc, :S_oc, :u_oc, :SIC, :h_ice, :T_ice_top, :T_sfc)
        d = maximum(abs, getproperty(s1, f) .- getproperty(s2, f))
        all_zero &= iszero(d)
        # Symbol → String so %-width formatting is unambiguous
        println(io, @sprintf("  max |Δ%-9s| = %g", String(f), d))
    end
    verdict = all_zero ? "IDENTICAL  (construction is deterministic)" :
                         "DIFFER     (global-state pollution or non-determinism)"
    println(io, "  → $verdict")
    flush(io)
    return all_zero
end

# Collect Test-1 extrema from already-materialized arrays / coupler fields.
function monitor_state_test1(ocean_model, ice_model, cs, Nz)
    T_oc      = interior_vec(ocean_model.tracers.T,       :, :, Nz)
    S_oc      = interior_vec(ocean_model.tracers.S,       :, :, Nz)
    u_oc      = interior_vec(ocean_model.velocities.u,    :, :, Nz)
    sic       = interior_vec(ice_model.ice_concentration, :, :, 1)
    h_ice     = interior_vec(ice_model.ice_thickness,     :, :, 1)
    T_ice_top = interior_vec(ice_model.ice_thermodynamics.top_surface_temperature, :, :, 1)
    T_sfc     = field_vec(cs.fields.T_sfc)
    F_sh      = field_vec(cs.fields.F_sh)
    F_lh      = field_vec(cs.fields.F_lh)
    return (; T_oc, S_oc, u_oc, sic, h_ice, T_ice_top, T_sfc, F_sh, F_lh)
end

function monitor_state_test2c(ocean_model, ice_model, cs, Nz)
    T_oc      = interior_vec(ocean_model.tracers.T, :, :, Nz)
    T_ice_top = interior_vec(ice_model.ice_thermodynamics.top_surface_temperature, :, :, 1)
    T_sfc     = field_vec(cs.fields.T_sfc)
    F_sh      = field_vec(cs.fields.F_sh)
    return (; T_oc, T_ice_top, T_sfc, F_sh)
end

function format_test1_line(ii, day, s)
    T_oc_lo, T_oc_hi           = finite_extrema(s.T_oc)
    S_oc_lo, S_oc_hi           = finite_extrema(s.S_oc)
    u_max                      = Float64(maximum(abs, s.u_oc))
    sic_lo, sic_hi             = finite_extrema(s.sic)
    h_lo, h_hi                 = finite_extrema(s.h_ice)
    T_ice_lo, T_ice_hi         = finite_extrema(s.T_ice_top)
    T_sfc_lo, T_sfc_hi         = finite_extrema(s.T_sfc)
    F_sh_lo, F_sh_hi           = finite_extrema(s.F_sh)
    F_lh_lo, F_lh_hi           = finite_extrema(s.F_lh)
    # Literal format only — do not concatenate format pieces with *.
    return @sprintf(
        "step %4d  day %6.2f | T_oc [%7.2f, %7.2f] C  S_oc [%5.2f, %5.2f] psu  u_oc_max %7.3f m/s | SIC [%.3f, %.3f]  h_ice [%.3f, %.3f] m | T_ice_top [%7.2f, %7.2f] C | T_sfc [%6.2f, %6.2f] K  F_sh [%8.2f, %8.2f] W/m2  F_lh [%8.2f, %8.2f] W/m2",
        ii, day,
        T_oc_lo, T_oc_hi,
        S_oc_lo, S_oc_hi,
        u_max,
        sic_lo, sic_hi,
        h_lo, h_hi,
        T_ice_lo, T_ice_hi,
        T_sfc_lo, T_sfc_hi,
        F_sh_lo, F_sh_hi,
        F_lh_lo, F_lh_hi,
    )
end

function format_test2c_line(ii, day, s)
    T_oc_lo, T_oc_hi     = finite_extrema(s.T_oc)
    T_ice_lo, T_ice_hi   = finite_extrema(s.T_ice_top)
    T_sfc_lo, T_sfc_hi   = finite_extrema(s.T_sfc)
    F_sh_lo, F_sh_hi     = finite_extrema(s.F_sh)
    return @sprintf(
        "2C step %4d  day %6.2f | T_oc [%7.2f, %7.2f] C | T_ice_top [%7.2f, %7.2f] C | T_sfc [%6.2f, %6.2f] K | F_sh [%8.2f, %8.2f] W/m2",
        ii, day,
        T_oc_lo, T_oc_hi,
        T_ice_lo, T_ice_hi,
        T_sfc_lo, T_sfc_hi,
        F_sh_lo, F_sh_hi,
    )
end

function monitor_has_nan_test1(s)
    return has_nan(s.T_oc) || has_nan(s.S_oc) || has_nan(s.u_oc) ||
           has_nan(s.sic) || has_nan(s.h_ice) || has_nan(s.T_ice_top) ||
           has_nan(s.T_sfc) || has_nan(s.F_sh) || has_nan(s.F_lh)
end

function monitor_has_nan_test2c(s)
    return has_nan(s.T_oc) || has_nan(s.T_ice_top) ||
           has_nan(s.T_sfc) || has_nan(s.F_sh)
end

# Smoke-test formatting before constructing CoupledSimulation (cheap catch for
# arg-count / literal-format mistakes that would otherwise burn a GPU node).
let
    dummy = fill(0.0, 4)
    s1 = (;
        T_oc = dummy, S_oc = dummy, u_oc = dummy, sic = dummy, h_ice = dummy,
        T_ice_top = dummy, T_sfc = dummy, F_sh = dummy, F_lh = dummy,
    )
    s2c = (; T_oc = dummy, T_ice_top = dummy, T_sfc = dummy, F_sh = dummy)
    line1 = format_test1_line(1, 0.0, s1)
    line2 = format_test2c_line(1, 0.0, s2c)
    @assert occursin("step", line1) && occursin("2C step", line2)
    @info "Format smoke-test OK"
end

# All runtime state lives inside `main` so assignments in `try`/`for` are ordinary
# locals (top-level soft-scope would otherwise leave nan_step/last_step stuck).
function main()
    # ─── Construct reference cs (baseline for Test 1 and Test 2B) ─────────────

    @info "──────────────────────────────────────────────────────────────────────"
    @info "Constructing reference CoupledSimulation..."
    cs_ref        = CoupledSimulation(CONFIG_FILE)
    snap_ref_init = snapshot_state(cs_ref)   # saved before any stepping

    # ─── Test 2A: back-to-back construction ───────────────────────────────────

    @info "══════ TEST 2A: back-to-back construction ══════"
    cs_2a   = CoupledSimulation(CONFIG_FILE)
    snap_2a = snapshot_state(cs_2a)

    open(joinpath(ARTIFACT_DIR, "test2A_comparison.txt"), "w") do io
        write_comparison(io,
            "Test 2A — back-to-back construction (expect all 0.0)",
            snap_ref_init, snap_2a)
    end
    write_comparison(stdout,
        "Test 2A — back-to-back construction (expect all 0.0)",
        snap_ref_init, snap_2a)

    # ─── Test 1: step-by-step extrema monitoring ──────────────────────────────

    @info "══════ TEST 1: step-by-step extrema (N_MONITOR=$N_MONITOR) ══════"

    nan_step  = nothing
    last_step = 0
    test1_log = joinpath(ARTIFACT_DIR, "test1_extrema.txt")
    log_io_1  = open(test1_log, "w")

    try
        ocean_model = cs_ref.model_sims.ocean_sim.ocean.model
        ice_model   = cs_ref.model_sims.ice_sim.ice.model
        Nz          = size(ocean_model.grid, 3)

        for ii in 1:N_MONITOR
            step!(cs_ref)
            last_step = ii
            day = Float64(cs_ref.t[]) / 86400

            s = monitor_state_test1(ocean_model, ice_model, cs_ref, Nz)
            line = format_test1_line(ii, day, s)
            println(log_io_1, line)
            flush(log_io_1)
            @info line

            if monitor_has_nan_test1(s)
                nan_step = ii
                @warn "NaN/Inf detected at step $ii (day $(round(day, digits=3))) — stopping Test 1"
                break
            end
        end
    catch e
        crashed_day = round(Float64(cs_ref.t[]) / 86400, digits = 3)
        @error "Test 1 threw at step $last_step (day $crashed_day)" exception = e
    finally
        close(log_io_1)
    end

    if isnothing(nan_step)
        @info "Test 1: completed $last_step steps with no NaN detected"
    else
        @info "Test 1: NaN first detected at step $nan_step  (day $(round(Float64(cs_ref.t[])/86400, digits=3)))"
    end

    # ─── Test 2B: fresh cs constructed after cs_ref has run or crashed ────────

    ref_day = round(Float64(cs_ref.t[]) / 86400, digits = 3)
    @info "══════ TEST 2B: fresh cs after Test 1 (cs_ref reached day $ref_day) ══════"
    cs_2b   = CoupledSimulation(CONFIG_FILE)
    snap_2b = snapshot_state(cs_2b)

    label_2b = "Test 2B — fresh cs after run (nonzero → global-state pollution; zero → same ICs, divergence is inside step!)"
    open(joinpath(ARTIFACT_DIR, "test2B_comparison.txt"), "w") do io
        println(io, "cs_ref ran to: day $ref_day")
        isnothing(nan_step) || println(io, "NaN first detected at step: $nan_step")
        write_comparison(io, label_2b, snap_ref_init, snap_2b)
    end
    write_comparison(stdout, label_2b, snap_ref_init, snap_2b)

    # ─── Test 2C: step cs_2b to see if it also crashes ───────────────────────
    #
    # Interpretation matrix:
    #   cs_ref crashes, cs_2b crashes at same step → crash is deterministic; if
    #     Test 2B ICs differ the IC difference is a red herring (crash is fundamental)
    #   cs_ref crashes, cs_2b crashes at different step → non-determinism inside step!
    #   cs_ref crashes, cs_2b stable → if Test 2B ICs differ: IC pollution explains
    #     divergence; if ICs same: non-determinism inside step!
    #   neither crashes → neither run is unstable at N_MONITOR steps (increase limit)

    @info "══════ TEST 2C: step cs_2b for N_MONITOR=$N_MONITOR steps ══════"

    nan_step_2c  = nothing
    last_step_2c = 0
    test2c_log = joinpath(ARTIFACT_DIR, "test2C_extrema.txt")
    log_io_2c  = open(test2c_log, "w")

    try
        ocean_model_2b = cs_2b.model_sims.ocean_sim.ocean.model
        ice_model_2b   = cs_2b.model_sims.ice_sim.ice.model
        Nz_2b          = size(ocean_model_2b.grid, 3)

        for ii in 1:N_MONITOR
            step!(cs_2b)
            last_step_2c = ii
            day = Float64(cs_2b.t[]) / 86400

            s = monitor_state_test2c(ocean_model_2b, ice_model_2b, cs_2b, Nz_2b)
            line = format_test2c_line(ii, day, s)
            println(log_io_2c, line)
            flush(log_io_2c)
            @info line

            if monitor_has_nan_test2c(s)
                nan_step_2c = ii
                @warn "Test 2C: NaN/Inf detected at step $ii (day $(round(day, digits=3))) — stopping"
                break
            end
        end
    catch e
        crashed_day_2c = round(Float64(cs_2b.t[]) / 86400, digits = 3)
        @error "Test 2C threw at step $last_step_2c (day $crashed_day_2c)" exception = e
    finally
        close(log_io_2c)
    end

    open(joinpath(ARTIFACT_DIR, "test2C_summary.txt"), "w") do io
        ref_result = isnothing(nan_step)    ? "stable for $last_step steps"  : "NaN at step $nan_step"
        new_result = isnothing(nan_step_2c) ? "stable for $last_step_2c steps" : "NaN at step $nan_step_2c"
        println(io, "Test 1  (cs_ref): $ref_result")
        println(io, "Test 2C (cs_2b):  $new_result")
        println(io)
        verdict = if isnothing(nan_step) && isnothing(nan_step_2c)
            "BOTH STABLE — crash not reproduced; increase N_MONITOR"
        elseif !isnothing(nan_step) && !isnothing(nan_step_2c) && nan_step == nan_step_2c
            "BOTH CRASH at same step — crash is deterministic and reproducible"
        elseif !isnothing(nan_step) && !isnothing(nan_step_2c)
            "BOTH CRASH but at different steps — non-determinism inside step! (same ICs, different crash steps)"
        elseif !isnothing(nan_step) && isnothing(nan_step_2c)
            "cs_ref crashes; cs_2b stable — if Test 2B ICs differ: IC pollution explains divergence; if ICs same: non-determinism inside step!"
        else
            "cs_ref stable; cs_2b crashes — unexpected ordering; check for measurement error"
        end
        println(io, "Interpretation: $verdict")
        flush(io)
    end

    if isnothing(nan_step_2c)
        @info "Test 2C: completed $last_step_2c steps with no NaN"
    else
        @info "Test 2C: NaN at step $nan_step_2c"
    end
    @info "Stability diagnostics complete. Artifacts written to $ARTIFACT_DIR"
    return nothing
end

main()
