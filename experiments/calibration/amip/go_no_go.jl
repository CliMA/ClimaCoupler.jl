# Iteration-1 go/no-go gate.
#
# WHY. A calibration's fate is largely decided by its FIRST ensemble, and that
# ensemble is cheap to read. The micro_edmf run spent six iterations and ~9
# hours to end up with residuals flat to within a percent, and both causes
# were visible at iteration 1:
#
#   - leverage ratios of 2.2-3.3 (residual well outside what the parameters
#     could move), and
#   - a reachable spread of 3.20 W/m^2 median against an 11.2 W/m^2 noise
#     floor, i.e. a parameter signal-to-noise of about 0.3.
#
# Run this as soon as iteration_001/G_ensemble.jld2 exists. It reports two
# quantities per observable and applies explicit thresholds.
#
#   LEVERAGE RATIO = residual / spread, both in noise-sigma units (from
#     campaign_reference/leverage.jl). Measured anchors on matched 5-degree
#     observations: below ~2 the observable improves freely; near 3 expect a
#     few percent; near 4 expect nothing, and reweighting cannot help because
#     the bias is out of the parameters' reach.
#
#   REACHABLE SPREAD = per-point standard deviation of G across members, in
#     PHYSICAL units. This is the question leverage's sigma-normalised spread
#     obscures: can these parameters move the field by an amount comparable to
#     the noise floor we told EKP to assume? If the median is well below the
#     floor, EKP has almost nothing to learn from and the run is a formality.
#
# ENSEMBLE CONTRACTION is also reported once a posterior exists (iteration 2
# onward). The micro_edmf run contracted 71x, the signature of an
# over-informative covariance.
#
# Usage, from the repository root:
#
#   julia --project=experiments/AMIP experiments/calibration/amip/go_no_go.jl \
#       <output_dir> [iteration]
#
# Environment overrides for the thresholds (defaults are the rlut_pigroups
# predictions): GNG_MAX_RATIO (3.0), GNG_MIN_SPREAD (5.0, W/m^2),
# GNG_MAX_CONTRACTION (10.0).
#
# Exits 0 on GO, 1 on NO-GO, so it can gate a submit script.

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2
import Statistics
import Printf: @printf, @sprintf

include(joinpath(@__DIR__, "campaign_reference", "leverage.jl"))

output_dir = length(ARGS) >= 1 ? ARGS[1] : error("Usage: go_no_go.jl <output_dir> [iteration]")
iteration = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

max_ratio = parse(Float64, get(ENV, "GNG_MAX_RATIO", "3.0"))
min_spread = parse(Float64, get(ENV, "GNG_MIN_SPREAD", "5.0"))
max_contraction = parse(Float64, get(ENV, "GNG_MAX_CONTRACTION", "10.0"))

println("="^76)
println("ITERATION-$iteration GO/NO-GO  -  $output_dir")
println("="^76)

failures = String[]
warnings = String[]

# --------------------------------------------------------------------------
# 1. Leverage, in noise-sigma units
# --------------------------------------------------------------------------
rows = leverage(output_dir; iteration)
println("\n1. LEVERAGE (noise-sigma units)")
@printf("   %-10s %8s %8s %8s %8s %7s\n",
        "observable", "n", "residual", "spread", "ratio", "reach")
for r in rows
    @printf("   %-10s %8d %8.3f %8.3f %8.2f %6.1f%%\n",
            r.name, r.n, r.residual, r.spread, r.ratio, 100 * r.reach)
    if !isfinite(r.ratio) || r.ratio > max_ratio
        push!(failures, "$(r.name): leverage ratio $(round(r.ratio; digits=2)) > $max_ratio")
    end
    if r.reach < 0.5
        push!(warnings,
              "$(r.name): observation inside the member envelope at only " *
              "$(round(100*r.reach; digits=1))% of points - the update extrapolates")
    end
end
@printf("   threshold: ratio <= %.1f\n", max_ratio)

# --------------------------------------------------------------------------
# 2. Reachable spread, in physical units
# --------------------------------------------------------------------------
ekp, _ = _leverage_load(output_dir)
y, sigma, _, blocks = _leverage_blocks(ekp, iteration)
g = JLD2.load_object(leverage_g_path(output_dir, iteration))

println("\n2. REACHABLE SPREAD (physical units - can the parameters move it?)")
@printf("   %-10s %10s %10s %10s %10s\n",
        "observable", "spread", "p90", "floor", "spread/floor")
for (name, rng) in blocks
    idx = [i for i in rng if i <= size(g, 1) && isfinite(y[i]) && all(isfinite, view(g, i, :))]
    isempty(idx) && continue
    per_point = [Statistics.std(view(g, i, :)) for i in idx]
    med = Statistics.median(per_point)
    p90 = Statistics.quantile(per_point, 0.9)
    floor_med = Statistics.median(sigma[idx])
    @printf("   %-10s %10.2f %10.2f %10.2f %10.2f\n",
            name, med, p90, floor_med, med / floor_med)
    if med < min_spread
        push!(failures,
              "$(name): reachable spread $(round(med; digits=2)) < $min_spread " *
              "- the parameters barely move this observable")
    end
end
@printf("   threshold: spread >= %.1f (physical units of the observable)\n", min_spread)

# --------------------------------------------------------------------------
# 3. Ensemble contraction, once a posterior exists
# --------------------------------------------------------------------------
println("\n3. ENSEMBLE CONTRACTION (parameter space, unconstrained)")
u = EKP.get_u(ekp)
if length(u) < 2
    println("   not available yet - needs at least one EKP update (iteration 2+)")
else
    s0 = Statistics.mean(vec(Statistics.std(u[1], dims = 2)))
    sf = Statistics.mean(vec(Statistics.std(u[end], dims = 2)))
    contraction = s0 / sf
    @printf("   prior %.4f -> current %.4f   contraction %.1fx (threshold <= %.0fx)\n",
            s0, sf, contraction, max_contraction)
    contraction > max_contraction && push!(failures,
        "ensemble contracted $(round(contraction; digits=1))x > $max_contraction " *
        "- the covariance is over-informative")
end

# --------------------------------------------------------------------------
# 4. Failed members
# --------------------------------------------------------------------------
n_members = size(g, 2)
dead = count(m -> !any(isfinite, view(g, :, m)), 1:n_members)
println("\n4. MEMBERS")
@printf("   %d of %d produced finite output\n", n_members - dead, n_members)
dead > 0 && push!(warnings, "$dead of $n_members members produced no finite output")

# --------------------------------------------------------------------------
println("\n", "="^76)
for w in warnings
    println("WARN   ", w)
end
if isempty(failures)
    println("VERDICT: GO")
    println("="^76)
    exit(0)
else
    for f in failures
        println("FAIL   ", f)
    end
    println("VERDICT: NO-GO  -  stop the run rather than spend the remaining iterations")
    println("="^76)
    exit(1)
end
