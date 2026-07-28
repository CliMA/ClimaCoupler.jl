# Generate a concise, plain-language report on a calibration run: its
# configuration and its results (parameters and covariance-weighted error across
# iterations, plus ensemble-spread convergence / collapse diagnosis).
#
# The report is (mostly) self-contained from the run's output_dir: parameter
# names and constrained transforms come from the saved `prior.jld2`, and the
# error/spread/parameter trajectories come from the latest `eki_file.jld2`.
# Passing the run's config .jl adds the target/covariance dates and iteration
# target (requires ClimaCoupler to be loadable, which it is in the AMIP env).
#
# Usage
# -----
#   # As a script (writes <output_dir>/calibration_report.md and prints to stdout):
#   julia --project=experiments/AMIP experiments/calibration/amip/calibration_report.jl \
#       <output_dir> [config_path] [label]
#
#   # As a function (after include):
#   calibration_report(output_dir; config_path = nothing, label = ..., outfile = ...)

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import ClimaCalibrate
import ClimaAnalysis
import JLD2
import Statistics
import LinearAlgebra
import YAML
import Printf: @sprintf

"""
    calibration_report(output_dir; config_path=nothing,
                       label=basename(rstrip(output_dir,'/')),
                       outfile=joinpath(output_dir,"calibration_report.md"))

Build a plain-language Markdown report for the calibration in `output_dir`,
print it, and (unless `isnothing(outfile)`) write it there. Returns the report
string.
"""
function calibration_report(
    output_dir;
    config_path = nothing,
    label = basename(rstrip(output_dir, '/')),
    outfile = joinpath(output_dir, "calibration_report.md"),
)
    isdir(output_dir) || error("output_dir does not exist: $output_dir")

    # --- locate iterations + load the latest ekp ---
    iters = sort(filter(d -> occursin(r"^iteration_\d+$", d), readdir(output_dir)))
    isempty(iters) && error("No iteration_* directories in $output_dir")
    ekp = JLD2.load_object(joinpath(output_dir, iters[end], "eki_file.jld2"))

    # --- prior (parameter names + constrained transform), from the saved file ---
    prior_path = joinpath(output_dir, iters[1], "prior.jld2")
    isfile(prior_path) || (prior_path = joinpath(output_dir, "prior.jld2"))
    prior = isfile(prior_path) ? JLD2.load_object(prior_path) : nothing
    pnames = isnothing(prior) ? nothing : PD.get_name(prior)

    # --- trajectories ---
    err = EKP.get_error(ekp)                     # one per completed update
    us = EKP.get_u(ekp)                          # unconstrained ensemble per iter
    spreads = [Statistics.mean(vec(Statistics.std(u, dims = 2))) for u in us]
    spread_ratio = spreads ./ spreads[1]
    n_param_iters = length(us)
    N_ens = EKP.get_N_ens(ekp)
    obs_dim = size(EKP.get_g(ekp)[1], 1)

    # constrained parameter ensemble mean per iteration (physical units)
    ϕmean = nothing
    if !isnothing(prior)
        ϕ = EKP.get_ϕ(prior, ekp)
        ϕmean = [vec(Statistics.mean(p, dims = 2)) for p in ϕ]
    end

    # --- observation variables + optional config context ---
    obs_names = nothing
    obs_path = joinpath(output_dir, "observation_vec.jld2")
    if isfile(obs_path)
        try
            o = JLD2.load_object(obs_path)
            obs_names = o[1].names
        catch
        end
    end

    cfg = nothing
    if !isnothing(config_path) && isfile(config_path)
        cfg = _load_config_context(config_path)
    end

    # observation analysis (data reduction + noise covariance); nothing if absent
    oa = _observation_analysis(output_dir)

    # --- collapse diagnosis ---
    collapsed = n_param_iters >= 2 && spread_ratio[end] < 0.01
    # a "healthy" contraction is gradual; flag a one-step crash specifically
    one_step_crash = n_param_iters >= 2 && spread_ratio[2] < 0.05

    # ===================== assemble the report =====================
    io = IOBuffer()
    p(args...) = println(io, args...)

    p("# Calibration report: ", label)
    p()
    p("- Output dir: `", output_dir, "`")
    p("- Iterations with parameters: ", n_param_iters,
      !isnothing(cfg) ? " (of $(cfg.n_iterations) targeted)" : "")
    p("- Ensemble members (N_ens): ", N_ens)
    p("- Observation dimension: ", obs_dim,
      !isnothing(obs_names) ? "  (variables: $(join(obs_names, ", ")))" : "")
    if !isnothing(cfg)
        p("- Target variables: ", join(cfg.short_names, ", "))
        p("- Calibration target date(s): ", _fmt_ranges(cfg.sample_ranges_unique))
        p("- Covariance estimated from: ", _fmt_ranges(cfg.covariance_ranges))
        p("- Spinup / averaging window (extend): ", cfg.spinup, " / ", cfg.extend)
    end
    p()

    # --- data reduction + covariance (if the observation vector is available) ---
    if !isnothing(oa)
        _reduction_section(io, oa, cfg)
        _covariance_section(io, oa)
    end

    # --- parameters ---
    p("## Parameters")
    p()
    if isnothing(prior)
        p("_prior.jld2 not found; parameter names unavailable._")
    else
        p("| parameter | initial | final | change |")
        p("|---|---|---|---|")
        for j in 1:length(pnames)
            i0 = ϕmean[1][j]
            i1 = ϕmean[end][j]
            pct = i0 == 0 ? "—" : @sprintf("%+.0f%%", 100 * (i1 - i0) / abs(i0))
            p("| `", pnames[j], "` | ", _sig(i0), " | ", _sig(i1), " | ", pct, " |")
        end
    end
    p()

    # --- error + spread trajectory ---
    p("## Error and ensemble spread across iterations")
    p()
    hdr = "| iter | cov-weighted error | spread (rel. to iter 1) |"
    sep = "|---|---|---|"
    if !isnothing(prior)
        hdr *= " " * join(["`$(n)` |" for n in pnames], " ")
        sep *= " " * join(["---|" for _ in pnames], " ")
    end
    p(hdr); p(sep)
    for i in 1:n_param_iters
        e = i <= length(err) ? @sprintf("%.4f", err[i]) : "—"
        row = "| $i | $e | $(@sprintf("%.3f", spread_ratio[i])) |"
        if !isnothing(prior)
            row *= " " * join([_sig(ϕmean[i][j]) * " |" for j in 1:length(pnames)], " ")
        end
        p(row)
    end
    p()

    # --- summary stats ---
    e1 = err[1]
    ebest = minimum(err)
    elast = err[end]
    p("## Summary")
    p()
    p("- Error: initial ", _sig(e1), ", final ", _sig(elast), ", best ", _sig(ebest),
      "  (best is ", @sprintf("%.0f%%", 100 * (e1 - ebest) / e1), " below initial).")
    p("- Ensemble spread: contracted to ",
      @sprintf("%.0f%%", 100 * spread_ratio[end]), " of its initial value.")
    verdict = collapsed ? (one_step_crash ?
        "**COLLAPSED** — spread crashed to <5% after a single update; the ensemble froze and later iterations only reflect model noise, not learning." :
        "**COLLAPSED** — spread fell below 1% of initial; parameters effectively frozen.") :
        "**No collapse** — spread contracted gradually, consistent with healthy convergence toward the posterior."
    p("- Convergence: ", verdict)
    p()

    # --- plain-language narrative ---
    p("## Plain-language summary")
    p()
    p(_narrative(label, cfg, n_param_iters, N_ens, obs_dim, obs_names,
                 pnames, ϕmean, err, spread_ratio, collapsed, one_step_crash))
    p()

    report = String(take!(io))
    print(report)
    if !isnothing(outfile)
        open(outfile, "w") do f
            write(f, report)
        end
        @info "Wrote report to $outfile"
    end
    return report
end

# ---- observation analysis (data reduction + noise covariance) ----

# Analyze the saved observation vector: reconstruct the per-variable fields (to
# report how the data was spatially reduced) and decompose the SVDplusD
# covariance diagonal (to report the effective noise). Returns `nothing` if no
# observation_vec.jld2 is present.
function _observation_analysis(output_dir)
    obs_path = joinpath(output_dir, "observation_vec.jld2")
    isfile(obs_path) || return nothing
    o = try
        JLD2.load_object(obs_path)
    catch
        return nothing
    end
    o1 = o[1]
    y = reduce(vcat, o1.samples)
    vars = ClimaCalibrate.ObservationRecipe.reconstruct_vars(o1)
    names = [ClimaAnalysis.short_name(v) for v in vars]
    finite = [count(isfinite, v.data) for v in vars]  # flatten order == vars order

    # SVDplusD diagonal: low-rank interannual part + diagonal floor (model error).
    cov = o1.covs[1]
    S = cov.svd_cov.S
    U = cov.svd_cov.U
    Vt = cov.svd_cov.Vt
    dlow = zeros(length(y))
    for k in 1:length(S)
        dlow .+= S[k] .* (U[:, k] .* Vt[k, :])
    end
    dfloor = cov.diag_cov.diag
    σtot = sqrt.(max.(dlow .+ dfloor, 0.0))
    σfloor = sqrt.(max.(dfloor, 0.0))

    # per-variable index ranges in the flattened vector
    offs = cumsum([0; finite])
    ok = offs[end] == length(y)
    pervar = ok ? [(name = names[i], v = vars[i], n = finite[i],
                    rng = (offs[i] + 1):offs[i + 1]) for i in 1:length(vars)] :
             nothing

    return (; o1, y, vars, names, finite, pervar, σtot, σfloor,
            rank = length(S), S)
end

function _reduction_section(io, oa, cfg)
    p(args...) = println(io, args...)
    p("## Data reduction for calibration")
    p()
    zonal = all(v -> ClimaAnalysis.has_latitude(v) && !ClimaAnalysis.has_longitude(v),
                oa.vars)
    if zonal
        p("- **Spatial: zonal (longitude) mean** — the longitude dimension is ",
          "averaged out (NaN-aware), so each field is reduced to latitude bands.")
    else
        p("- Spatial: fields kept on their lon/lat grid (no zonal averaging detected).")
    end
    p()
    p("| variable | reduced grid | calibration points |")
    p("|---|---|---|")
    for pv in (isnothing(oa.pervar) ? [] : oa.pervar)
        v = pv.v
        dimdesc = String[]
        if ClimaAnalysis.has_altitude(v)
            alts = ClimaAnalysis.altitudes(v)
            altstr = join(Int.(round.(alts)), "/")
            push!(dimdesc, "$(length(alts)) levels ($(altstr) m)")
        end
        if ClimaAnalysis.has_latitude(v)
            push!(dimdesc, "$(length(ClimaAnalysis.latitudes(v))) latitudes")
        end
        p("| `", pv.name, "` | ", join(dimdesc, " × "), " | ", pv.n, " |")
    end
    p("| **total** | | **", length(oa.y), "** |")
    p()
    if zonal && !isnothing(cfg) && !isnothing(cfg.nlon)
        nlat = length(ClimaAnalysis.latitudes(oa.vars[1]))
        full = cfg.nlon * nlat * sum(v -> ClimaAnalysis.has_altitude(v) ?
                                     length(ClimaAnalysis.altitudes(v)) : 1, oa.vars)
        p("- The model grid has ~", cfg.nlon, " longitudes per latitude, so the ",
          "zonal mean collapses each latitude's ", cfg.nlon,
          " values to one. Combined with level selection this is ",
          length(oa.y), " points — roughly **",
          @sprintf("%d×", round(Int, full / length(oa.y))),
          " fewer** than the full lon×lat grid (~", full, ").")
        p("- Why it matters: the SVDplusD covariance only softens a handful of ",
          "interannual modes and otherwise treats each point as independent; ",
          "cutting spatially-correlated points to the true large-scale degrees of ",
          "freedom is what keeps the inverse from over-informing and collapsing the ensemble.")
    end
    p()
end

function _covariance_section(io, oa)
    p(args...) = println(io, args...)
    p("## Observation covariance (noise model)")
    p()
    p("- Type: **SVDplusD** — a low-rank interannual covariance (estimated from the ",
      "spread across the covariance dates) plus a diagonal structural-error floor.")
    p("- SVD rank: ", oa.rank, "  (singular values: ",
      join([_sig(s) for s in oa.S], ", "), ")")
    p()
    p("Effective noise σ = √diag(Γ), per variable (median over points):")
    p()
    p("| variable | data \\|y\\| | noise σ | σ as % of field | floor σ (model error) |")
    p("|---|---|---|---|---|")
    for pv in (isnothing(oa.pervar) ? [] : oa.pervar)
        r = pv.rng
        yv = Statistics.median(abs.(oa.y[r]))
        st = Statistics.median(oa.σtot[r])
        sf = Statistics.median(oa.σfloor[r])
        pct = yv == 0 ? "—" : @sprintf("%.0f%%", 100 * st / yv)
        pctf = yv == 0 ? "—" : @sprintf("%.0f%%", 100 * sf / yv)
        p("| `", pv.name, "` | ", _sig(yv), " | ", _sig(st), " | ", pct,
          " | ", _sig(sf), " (", pctf, ") |")
    end
    p()
    p("- Interpretation: where the noise σ is comparable to or larger than the ",
      "model–obs mismatch, that variable carries little usable signal (EKP treats ",
      "the mismatch as noise); where σ is well below the mismatch, it drives the fit. ",
      "The diagonal floor encodes the structural-error allowance a perfect ",
      "parameter set is expected to leave.")
    p()
end

# ---- helpers ----

_sig(x) = abs(x) >= 100 ? @sprintf("%.0f", x) :
          abs(x) >= 1 ? @sprintf("%.3g", x) : @sprintf("%.3g", x)

function _fmt_ranges(ranges)
    isempty(ranges) && return "—"
    fmt(d) = string(d)[1:10]  # yyyy-mm-dd
    parts = [first(r) == last(r) ? fmt(first(r)) : "$(fmt(first(r)))..$(fmt(last(r)))"
             for r in ranges]
    return join(parts, ", ")
end

# Load just the context we want from a config .jl without launching anything.
function _load_config_context(config_path)
    m = Module(:CalibReportCfg)
    Core.eval(m, :(import ClimaCoupler, ClimaCoupler.CalibrationTools))
    Core.eval(m, :(import EnsembleKalmanProcesses as EKP))
    Core.eval(m, :(import EnsembleKalmanProcesses.ParameterDistributions as PD))
    Core.eval(m, :(import Dates))
    Base.include(m, abspath(config_path))
    C = Core.eval(m, :CALIBRATE_CONFIG)
    covr = Core.eval(m, :COVARIANCE_DATE_RANGES)
    sr = C.sample_date_ranges
    # model longitude count (for the zonal-mean reduction factor), from h_elem
    nlon = nothing
    try
        cd = YAML.load_file(C.config_file)
        haskey(cd, "h_elem") && (nlon = cd["h_elem"] * 4 * 3)
    catch
    end
    return (
        short_names = C.short_names,
        n_iterations = C.n_iterations,
        sample_ranges_unique = unique(sr),
        covariance_ranges = covr,
        spinup = C.spinup,
        extend = C.extend,
        nlon = nlon,
    )
end

function _narrative(label, cfg, niter, nens, obsdim, obsnames, pnames, ϕmean,
                    err, spread_ratio, collapsed, one_step_crash)
    io = IOBuffer()
    vars = !isnothing(obsnames) ? join(obsnames, "+") :
           (!isnothing(cfg) ? join(cfg.short_names, "+") : "the target field(s)")
    target = !isnothing(cfg) ? _fmt_ranges(cfg.sample_ranges_unique) : "a fixed target"
    print(io, "This calibration fit $(nens) ensemble members to $(obsdim) data points ",
          "of $(vars) (target: $(target)) over $(niter) iterations")
    !isnothing(cfg) && niter < cfg.n_iterations &&
        print(io, " (short of the $(cfg.n_iterations) targeted, typically because the workers hit their walltime)")
    print(io, ". ")

    e1, elast, ebest = err[1], err[end], minimum(err)
    redpct = round(100 * (e1 - ebest) / e1)
    if collapsed
        print(io, "The ensemble spread ",
              one_step_crash ? "crashed after the first update" : "fell below 1% of its initial value",
              ", so the parameters froze and the error (initial $(round(e1,sigdigits=3)), ",
              "wandering around $(round(Statistics.mean(err),sigdigits=3))) reflects model/weather ",
              "noise rather than learning. This calibration did not meaningfully improve the fit. ")
    else
        trend = elast < e1 ? "decreased" : "stayed roughly level"
        print(io, "The ensemble spread contracted gradually (to $(round(100*spread_ratio[end]))% of ",
              "initial), i.e. it converged rather than collapsing. The covariance-weighted error ",
              "$(trend) from $(round(e1,sigdigits=3)) to $(round(elast,sigdigits=3)) ",
              "(best $(round(ebest,sigdigits=3)), $(redpct)% below the start). ")
    end

    if !isnothing(pnames) && !collapsed
        moves = String[]
        for j in 1:length(pnames)
            i0, i1 = ϕmean[1][j], ϕmean[end][j]
            push!(moves, "$(pnames[j]) $(_sig(i0)) → $(_sig(i1))")
        end
        print(io, "Recovered parameters: ", join(moves, "; "), ".")
    end
    return String(take!(io))
end

# ---- CLI entrypoint ----
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: calibration_report.jl <output_dir> [config_path] [label]")
    output_dir = ARGS[1]
    config_path = length(ARGS) >= 2 && !isempty(ARGS[2]) ? ARGS[2] : nothing
    label = length(ARGS) >= 3 ? ARGS[3] : basename(rstrip(output_dir, '/'))
    calibration_report(output_dir; config_path, label)
end
