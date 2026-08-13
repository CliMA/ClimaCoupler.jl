# Per-observable leverage diagnostic: can this parameter set move this
# observable at all?
#
# For one iteration it reports, per observable, in units of the assumed
# noise sigma:
#   residual  RMS of (mean(G) - y) / sigma, how far the ensemble mean sits
#             from the observation.
#   spread    mean over members of RMS of (G_member - mean(G)) / sigma, how
#             much the observable moves when the parameters span their
#             ensemble. This is the leverage.
#   ratio     residual / spread. It is the excursion, in units of the
#             current ensemble width, needed to close the residual. A ratio
#             far above 1 means the bias is out of reach of these
#             parameters, and no reweighting of the observable will help.
#   reach     fraction of points where y lies inside the member envelope.
#             A low value means the update extrapolates.
#
# Measured anchor points for the ratio, all on the same October 2010
# lwp/clt/swcre observations at 5 degrees:
#   1.5-1.8  lwp. Both runs improved it.
#   2.6      swcre. Both runs improved it, by 0.09 sigma.
#   2.7      clt in relax, with the cloud-fraction shape free. It improved
#            by 0.066 sigma, that is, about 10 percent.
#   3.7      clt in vterm and clt_only_vterm, with the shape frozen. It
#            DEGRADED, 0.620 to 0.654 over six iterations.
# So the ratio predicts difficulty, not impossibility. Below about 2 the
# observable improves freely. Near 3 expect a few percent. Near 4 expect
# nothing, and reweighting will not help because the residual is out of
# reach of the parameters; a new parameter is the only fix.
#
# Compare spread across runs only at matched priors. A wider prior widens
# the G spread by itself, so relax's larger clt spread (0.241 against
# vterm's 0.166) mixes its three extra parameters with its wider q_liq and
# timescale priors.
#
# Usage, from the repository root:
#
#   julia --project=experiments/AMIP experiments/calibration/amip/leverage.jl \
#       <output_dir> [iteration]
#
#   # or, after include:
#   rows = leverage(output_dir; iteration = 1)
#   print_leverage(rows; label = "clt_only")

import ClimaAnalysis
import ClimaCalibrate
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses as LeverageEKP
import EnsembleKalmanProcesses.ParameterDistributions as LeverageEKPPD
import JLD2
import LinearAlgebra
import Statistics
import Printf: @printf, @sprintf

"Iterations that have a G ensemble, plus the ekp holding the full history."
function _leverage_load(output_dir)
    iters = sort(filter(d -> occursin(r"^iteration_\d+$", d), readdir(output_dir)))
    isempty(iters) && error("No iteration_* directories in $output_dir")
    ekp_path = nothing
    for it in reverse(iters)
        p = joinpath(output_dir, it, "eki_file.jld2")
        isfile(p) && (ekp_path = p; break)
    end
    isnothing(ekp_path) && error("No eki_file.jld2 in any iteration of $output_dir")
    g_iters = Int[]
    for it in iters
        n = parse(Int, match(r"\d+$", it).match)
        isfile(joinpath(output_dir, it, "G_ensemble.jld2")) && push!(g_iters, n)
    end
    isempty(g_iters) && error("No G_ensemble.jld2 in any iteration of $output_dir")
    return JLD2.load_object(ekp_path), sort(g_iters)
end

"Path to an iteration's G ensemble."
leverage_g_path(output_dir, i) =
    joinpath(output_dir, @sprintf("iteration_%03d", i), "G_ensemble.jld2")

"""
Observation vector, noise sigma, and one (name, index range) pair per
observable, for the minibatch that `iteration` was graded against. The
ranges follow the flattened order that the observation map writes.
"""
function _leverage_blocks(ekp, iteration)
    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs = ClimaCalibrate.get_observations_for_nth_iteration(obs_series, iteration)
    y = mapreduce(EKP.get_obs, vcat, minibatch_obs)
    sigma = mapreduce(vcat, minibatch_obs) do obs
        cov = EKP.get_obs_noise_cov(obs)
        [sqrt(abs(cov[i, i])) for i in 1:size(cov, 1)]
    end
    covs = [EKP.get_obs_noise_cov(obs) for obs in minibatch_obs]
    vars = mapreduce(ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, minibatch_obs)
    blocks = Tuple{String, UnitRange{Int}}[]
    offset = 0
    for v in vars
        len = length(ClimaAnalysis.flatten(v).data)
        push!(blocks, (ClimaAnalysis.short_name(v), (offset + 1):(offset + len)))
        offset += len
    end
    return y, sigma, covs, blocks
end

"""
    leverage(output_dir; iteration = nothing, whiten = false)

Per-observable residual, ensemble spread, their ratio, and reachability at
`iteration` (default: the first iteration with a G ensemble, which measures
the prior ensemble and so predicts the run before it spends its budget).

Set `whiten = true` to also report the residual and spread under the full
covariance block instead of its diagonal. That is the quantity the update
actually sees, and it differs a lot when the noise floor is spatially
correlated: a smooth bias is cheap under a long decorrelation length.
Whitening factorizes a dense block per observable, so it costs memory.

Returns a vector of NamedTuples, one per observable.
"""
function leverage(output_dir; iteration = nothing, whiten = false)
    ekp, g_iters = _leverage_load(output_dir)
    it = isnothing(iteration) ? first(g_iters) : iteration
    it in g_iters || error("Iteration $it has no G_ensemble.jld2 (available: $g_iters)")
    g = JLD2.load_object(
        joinpath(output_dir, @sprintf("iteration_%03d", it), "G_ensemble.jld2"),
    )
    y, sigma, covs, blocks = _leverage_blocks(ekp, it)
    # One covariance per minibatch member; whitening needs a single block.
    whiten &&
        length(covs) != 1 &&
        error("whiten = true needs a single observation in the minibatch")

    rows = []
    for (name, rng) in blocks
        idx = [
            i for i in rng if i <= size(g, 1) &&
                sigma[i] > 0 &&
                isfinite(y[i]) &&
                all(isfinite, view(g, i, :))
        ]
        isempty(idx) && continue
        gm = vec(Statistics.mean(view(g, idx, :), dims = 2))
        s = sigma[idx]
        residual = sqrt(Statistics.mean(abs2, (gm .- y[idx]) ./ s))
        per_member = [
            sqrt(Statistics.mean(abs2, (view(g, idx, m) .- gm) ./ s)) for m in 1:size(g, 2)
        ]
        spread = Statistics.mean(per_member)
        lo = vec(minimum(view(g, idx, :), dims = 2))
        hi = vec(maximum(view(g, idx, :), dims = 2))
        # lo and hi are local to this block; only y and g take global indices.
        reach = count(k -> lo[k] <= y[idx[k]] <= hi[k], eachindex(idx)) / length(idx)

        wres = wspread = NaN
        if whiten
            block = Matrix(covs[1][idx, idx])
            F = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(block))
            w(v) = sqrt(LinearAlgebra.dot(v, F \ v) / length(v))
            wres = w(gm .- y[idx])
            wspread = Statistics.mean(w(view(g, idx, m) .- gm) for m in 1:size(g, 2))
        end

        push!(
            rows,
            (;
                name,
                iteration = it,
                n = length(idx),
                residual,
                spread,
                ratio = spread > 0 ? residual / spread : Inf,
                reach,
                spread_min = minimum(per_member),
                spread_max = maximum(per_member),
                whitened_residual = wres,
                whitened_spread = wspread,
            ),
        )
    end
    return rows
end

"""
    weather_floor_valid(ekp)

Whether the spread-contraction estimate of the weather floor can be trusted.

That estimate treats the FINAL ensemble spread in G as pure weather, which
holds only if the parameter spread contracted enough that the parameter
contribution is negligible. Require contraction below 0.3.

The so_jan run contracted to 0.9926 and the estimate returned a signal of
about zero by construction, reporting an S/N of 0.70 that was an artefact.
"""
function weather_floor_valid(ekp; threshold = 0.3)
    us = LeverageEKP.get_u(ekp)
    length(us) < 2 && return (false, NaN)
    u1, uN = first(us), last(us)
    f = Statistics.mean([
        Statistics.std(uN[k, :]) / Statistics.std(u1[k, :]) for k in 1:size(u1, 1)
    ])
    return (f < threshold, f)
end

"Print the table `leverage` returns, one line per observable."
function print_leverage(rows; label = "", io = stdout)
    isempty(rows) && (println(io, "no observables with usable points"); return)
    isempty(label) || println(io, label, "  (iteration ", first(rows).iteration, ")")
    has_w = any(r -> isfinite(r.whitened_residual), rows)
    @printf(
        io,
        "  %-8s %7s %8s %8s %7s %7s%s\n",
        "obs",
        "points",
        "residual",
        "spread",
        "ratio",
        "reach",
        has_w ? "   wh.resid  wh.spread" : ""
    )
    for r in rows
        @printf(
            io,
            "  %-8s %7d %8.3f %8.3f %7.1f %7.2f",
            r.name,
            r.n,
            r.residual,
            r.spread,
            r.ratio,
            r.reach
        )
        has_w && @printf(io, " %10.3f %10.3f", r.whitened_residual, r.whitened_spread)
        println(io)
    end
    for r in rows
        r.spread < 0.2 && println(
            io,
            "  NOTE ",
            r.name,
            ": spread ",
            @sprintf("%.3f", r.spread),
            " sigma, the parameters barely move it",
        )
        r.ratio > 3.5 && println(
            io,
            "  NOTE ",
            r.name,
            ": ratio ",
            @sprintf("%.1f", r.ratio),
            ", out of reach of these parameters. Expect no improvement, and ",
            "reweighting will not help. Add a parameter this observable ",
            "responds to.",
        )
        2.5 < r.ratio <= 3.5 && println(
            io,
            "  NOTE ",
            r.name,
            ": ratio ",
            @sprintf("%.1f", r.ratio),
            ", hard to reach. Expect a few percent at best.",
        )
    end
    return nothing
end

"""
    leverage_by_parameter(output_dir)

Attribute iteration 1's ensemble spread to individual parameters, at no
extra cost.

TransformUnscented builds its 2p+1 members as the centre plus plus/minus
perturbations along the columns of the Cholesky factor of the prior
covariance. Independent priors make that factor diagonal, so each member
perturbs exactly ONE parameter and iteration 1 is already a one-at-a-time
sweep. This function finds which parameter each member moved, then reports
that member's RMS change in each observable, in noise units.

Only iteration 1 works. After the first update the parameter covariance is
no longer diagonal and the sigma points mix parameters. The function
errors rather than guess if the members are not one-at-a-time.

Returns (rows, unattributed) where rows have fields `parameter`,
`observable`, `delta` (RMS change in noise units) and `n_members`, and
`unattributed` is the quadrature check: sqrt of the summed squares of the
per-parameter deltas should approach the aggregate spread from
`leverage` when the response is locally linear.
"""
function leverage_by_parameter(output_dir)
    ekp, g_iters = _leverage_load(output_dir)
    first(g_iters) == 1 ||
        error("leverage_by_parameter needs iteration 1, found $(g_iters)")
    prior_path = joinpath(output_dir, "iteration_001", "prior.jld2")
    isfile(prior_path) || (prior_path = joinpath(output_dir, "prior.jld2"))
    isfile(prior_path) || error("prior.jld2 not found in $output_dir")
    prior = JLD2.load_object(prior_path)
    pnames = LeverageEKPPD.get_name(prior)
    u = EKP.get_u(ekp, 1)                       # unconstrained, p x N_ens
    g = JLD2.load_object(joinpath(output_dir, "iteration_001", "G_ensemble.jld2"))
    y, sigma, _, blocks = _leverage_blocks(ekp, 1)

    # The centre member differs from no other member in only one coordinate.
    umean = vec(Statistics.mean(u, dims = 2))
    dist = [maximum(abs.(u[:, m] .- umean)) for m in 1:size(u, 2)]
    centre = argmin(dist)
    moved = Dict{Int, Vector{Int}}()             # parameter index -> members
    for m in 1:size(u, 2)
        m == centre && continue
        d = abs.(u[:, m] .- u[:, centre])
        scale = maximum(d)
        scale == 0 && continue
        big = findall(>(0.01 * scale), d)
        length(big) == 1 || error(
            "member $m moved $(length(big)) parameters at once, so the " *
            "sigma points are not one-at-a-time and attribution is invalid",
        )
        push!(get!(moved, big[1], Int[]), m)
    end

    gc = view(g, :, centre)
    rows = []
    for (k, members) in sort(collect(moved))
        for (name, rng) in blocks
            idx = [
                i for i in rng if i <= size(g, 1) &&
                    sigma[i] > 0 &&
                    isfinite(y[i]) &&
                    isfinite(gc[i]) &&
                    all(m -> isfinite(g[i, m]), members)
            ]
            isempty(idx) && continue
            delta = Statistics.mean(
                sqrt(Statistics.mean(abs2, (view(g, idx, m) .- gc[idx]) ./ sigma[idx]))
                for m in members
            )
            push!(
                rows,
                (;
                    parameter = pnames[k],
                    observable = name,
                    delta,
                    n_members = length(members),
                ),
            )
        end
    end
    return rows
end

"Print `leverage_by_parameter` output as an observable-by-parameter table."
function print_leverage_by_parameter(rows; io = stdout)
    isempty(rows) && (println(io, "nothing to attribute"); return)
    obs = unique(r.observable for r in rows)
    params = unique(r.parameter for r in rows)
    @printf(io, "  per-parameter response, RMS change in noise units\n")
    @printf(io, "  %-58s", "parameter")
    for o in obs
        @printf(io, "%10s", o)
    end
    println(io)
    for p in params
        @printf(io, "  %-58s", p)
        for o in obs
            k = findfirst(r -> r.parameter == p && r.observable == o, rows)
            isnothing(k) ? @printf(io, "%10s", "-") : @printf(io, "%10.3f", rows[k].delta)
        end
        println(io)
    end
    @printf(io, "  %-58s", "quadrature sum (compare to aggregate spread)")
    for o in obs
        s = sqrt(sum(r.delta^2 for r in rows if r.observable == o))
        @printf(io, "%10.3f", s)
    end
    println(io)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: leverage.jl <output_dir> [iteration]")
    output_dir = ARGS[1]
    iteration = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nothing
    rows = leverage(output_dir; iteration)
    print_leverage(rows; label = basename(rstrip(output_dir, '/')))
    if isnothing(iteration) || iteration == 1
        println()
        try
            print_leverage_by_parameter(leverage_by_parameter(output_dir))
        catch e
            @warn "Per-parameter attribution unavailable" exception = e
        end
    end
end
