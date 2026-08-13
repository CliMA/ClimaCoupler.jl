# PDF summary report for a calibration run, built entirely from the run's
# output_dir (eki_file.jld2, prior.jld2, G_ensemble.jld2, observation_vec.jld2,
# and the PNGs the calibration already produces).
#
# Pages:
#   1. Parameters and error across iterations (the run's own
#      constrained_params_and_error.png).
#   2. Bias maps at the first and last iteration, side by side.
#   3. Per-month residual RMS across iterations (lead-time drift check).
#
# Usage, from the repository root:
#
#   julia --project=experiments/AMIP experiments/calibration/amip/calibration_report_pdf.jl \
#       <output_dir> [label]

import CairoMakie
import ClimaAnalysis
import ClimaCalibrate
import EnsembleKalmanProcesses as EKP
import JLD2
import PNGFiles
import Poppler_jll
import Statistics
import Printf: @sprintf

const PAGE_W = 1400

"Iterations that have a G ensemble, and the ekp with the full history."
function _load_run(output_dir)
    iters = sort(filter(d -> occursin(r"^iteration_\d+$", d), readdir(output_dir)))
    isempty(iters) && error("No iteration_* directories in $output_dir")
    ekp_path = nothing
    for it in reverse(iters)
        p = joinpath(output_dir, it, "eki_file.jld2")
        isfile(p) && (ekp_path = p; break)
    end
    isnothing(ekp_path) && error("No eki_file.jld2 in any iteration of $output_dir")
    ekp = JLD2.load_object(ekp_path)
    prior_path = joinpath(output_dir, first(iters), "prior.jld2")
    isfile(prior_path) || (prior_path = joinpath(output_dir, "prior.jld2"))
    isfile(prior_path) || error("prior.jld2 not found in $output_dir")
    prior = JLD2.load_object(prior_path)
    g_iters = Int[]
    for it in iters
        n = parse(Int, match(r"\d+$", it).match)
        isfile(joinpath(output_dir, it, "G_ensemble.jld2")) && push!(g_iters, n)
    end
    return ekp, prior, sort(g_iters)
end

_g_path(output_dir, i) =
    joinpath(output_dir, @sprintf("iteration_%03d", i), "G_ensemble.jld2")

"Row means of the G ensemble. NaN entries from failed members are skipped."
function _g_mean(g_ensemble)
    return map(eachrow(g_ensemble)) do row
        finite = filter(isfinite, row)
        isempty(finite) ? NaN : Statistics.mean(finite)
    end
end

"""
Per-observable, per-month residual RMS in noise units. Months are identified
by tagging a copy of each reconstructed variable with its time index and
flattening it, so the grouping is exact under any dimension order and NaN mask.
"""
function _obs_rms_by_month(ekp, g_ensemble, iteration)
    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs = ClimaCalibrate.get_observations_for_nth_iteration(obs_series, iteration)
    y = mapreduce(EKP.get_obs, vcat, minibatch_obs)
    sigma = mapreduce(vcat, minibatch_obs) do obs
        cov = EKP.get_obs_noise_cov(obs)
        [sqrt(abs(cov[i, i])) for i in 1:size(cov, 1)]
    end
    vars = mapreduce(ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, minibatch_obs)
    gm = _g_mean(g_ensemble)
    n = min(length(gm), length(y), length(sigma))
    out = []
    offset = 0
    for v in vars
        name = ClimaAnalysis.short_name(v)
        dimnames = collect(keys(v.dims))
        tpos = findfirst(
            d -> ClimaAnalysis.conventional_dim_name(d) in ("time", "date"),
            dimnames,
        )
        isnothing(tpos) && continue
        nt = length(v.dims[dimnames[tpos]])
        tags = Array{Float64}(undef, size(v.data))
        for (t, sl) in enumerate(eachslice(tags; dims = tpos))
            sl .= t
        end
        tags[isnan.(v.data)] .= NaN
        vtag = ClimaAnalysis.OutputVar(v.attributes, v.dims, v.dim_attributes, tags)
        flat_tags = ClimaAnalysis.flatten(vtag).data
        len = length(ClimaAnalysis.flatten(v).data)
        idx = (offset + 1):(offset + len)
        offset += len
        rms = map(1:nt) do t
            r = [
                (gm[i] - y[i]) / sigma[i] for (k, i) in enumerate(idx) if
                i <= n && flat_tags[k] == t && isfinite(gm[i]) && sigma[i] > 0
            ]
            isempty(r) ? NaN : sqrt(Statistics.mean(abs2, r))
        end
        push!(out, (; name, rms))
    end
    return out
end

# ---------------------------------------------------------------- pages

"""
The same per-panel EKP.Visualize plots that make the run's wide
constrained_params_and_error.png strip, reflowed into a grid so the page
matches the report's other pages.
"""
function _page_params_error(ekp, prior, label; ncols = 4)
    dim_size = sum(length.(EKP.batch(prior)))
    panels = dim_size + 2
    nrows = cld(panels, ncols)
    fig = CairoMakie.Figure(size = (PAGE_W, 360 * nrows + 50))
    CairoMakie.Label(
        fig[0, 1:ncols],
        "$label — parameters and error";
        fontsize = 24,
        font = :bold,
        tellwidth = false,
    )
    grid(k) = (fld(k - 1, ncols) + 1, mod1(k, ncols))
    for i in 1:dim_size
        r, c = grid(i)
        EKP.Visualize.plot_ϕ_over_iters(fig[r, c], ekp, prior, i)
    end
    r, c = grid(dim_size + 1)
    EKP.Visualize.plot_error_over_iters(fig[r, c], ekp, error_metric = "loss")
    r, c = grid(dim_size + 2)
    EKP.Visualize.plot_error_over_time(fig[r, c], ekp, error_metric = "loss")
    return fig
end

"First and last bias maps side by side for direct comparison."
function _page_bias_pair(paths_titles, label)
    imgs = [PNGFiles.load(p) for (p, _) in paths_titles]
    ncol = length(imgs)
    h, w = size(imgs[1])
    colw = PAGE_W / ncol
    fig = CairoMakie.Figure(
        size = (PAGE_W, round(Int, h * colw / w) + 40),
        figure_padding = 0,
    )
    CairoMakie.Label(
        fig[0, 1:ncol],
        "$label — bias maps";
        fontsize = 26,
        font = :bold,
        tellwidth = false,
    )
    for (c, (img, (_, title))) in enumerate(zip(imgs, paths_titles))
        ax = CairoMakie.Axis(fig[1, c])
        CairoMakie.hidedecorations!(ax)
        CairoMakie.hidespines!(ax)
        CairoMakie.image!(ax, CairoMakie.rotr90(img))
        CairoMakie.tightlimits!(ax)
        CairoMakie.text!(
            ax,
            0.5,
            0.99;
            text = title,
            space = :relative,
            align = (:center, :top),
            fontsize = 20,
            font = :bold,
        )
    end
    CairoMakie.colgap!(fig.layout, 8)
    return fig
end

function _page_monthly(ekp, output_dir, g_iters, label)
    rows = Dict{String, Dict{Int, Vector{Float64}}}()
    nt_max = 0
    for i in g_iters
        g = JLD2.load_object(_g_path(output_dir, i))
        for e in _obs_rms_by_month(ekp, g, i)
            get!(rows, e.name, Dict{Int, Vector{Float64}}())[i] = e.rms
            nt_max = max(nt_max, length(e.rms))
        end
    end
    nt_max >= 1 || return nothing
    names = sort(collect(keys(rows)))
    title =
        (
            nt_max == 1 ? "$label — per-variable residual" :
            "$label — per-month residual (month 1 = first sample month)"
        ) * ", RMS of (G - y)/noise: 1 = fit to the assumed noise"
    fig = CairoMakie.Figure(size = (PAGE_W, 380))
    CairoMakie.Label(
        fig[0, 1:length(names)],
        title;
        fontsize = 22,
        font = :bold,
        tellwidth = false,
    )
    for (c, name) in enumerate(names)
        ax = CairoMakie.Axis(
            fig[1, c],
            title = name,
            xlabel = "iteration",
            ylabel = "RMS [x noise]",
        )
        for t in 1:nt_max
            pts = sort([
                (i, rms[t]) for
                (i, rms) in rows[name] if t <= length(rms) && isfinite(rms[t])
            ])
            isempty(pts) && continue
            CairoMakie.scatterlines!(ax, first.(pts), last.(pts); label = "month $t")
        end
        c == 1 && nt_max > 1 && CairoMakie.axislegend(ax; labelsize = 10)
    end
    return fig
end

# ---------------------------------------------------------------- assembly

"""
    calibration_report_pdf(output_dir; label, outfile)

Build the PDF report for the run in `output_dir` and write it to `outfile`
(default `<output_dir>/calibration_report.pdf`). Pages that cannot be built
(missing inputs) are skipped with a warning. Returns `outfile`.
"""
function calibration_report_pdf(
    output_dir;
    label = basename(rstrip(output_dir, '/')),
    outfile = joinpath(output_dir, "calibration_report.pdf"),
)
    ekp, prior, g_iters = _load_run(output_dir)
    isempty(g_iters) && error("No G_ensemble.jld2 in any iteration of $output_dir")

    tmp = mktempdir()
    pages = String[]
    function addpage(name, build)
        fig = try
            build()
        catch e
            @warn "Skipping report page" name exception = e
            return
        end
        isnothing(fig) && return
        path = joinpath(tmp, name * ".pdf")
        CairoMakie.save(path, fig)
        push!(pages, path)
    end

    addpage("1_params_error", () -> _page_params_error(ekp, prior, label))

    bias = [
        (
            joinpath(output_dir, @sprintf("iteration_%03d", i), "bias_sample_dates.png"),
            "iteration $i",
        ) for i in unique([first(g_iters), last(g_iters)])
    ]
    filter!(pt -> isfile(pt[1]), bias)
    isempty(bias) || addpage("2_bias_first_last", () -> _page_bias_pair(bias, label))

    addpage("3_monthly", () -> _page_monthly(ekp, output_dir, g_iters, label))

    isempty(pages) && error("No report pages could be built for $output_dir")
    run(`$(Poppler_jll.pdfunite()) $(pages) $(outfile)`)
    @info "Wrote report" outfile n_pages = length(pages)
    return outfile
end

# ---------------------------------------------------------------- CLI

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: calibration_report_pdf.jl <output_dir> [label]")
    output_dir = ARGS[1]
    label = length(ARGS) >= 2 ? ARGS[2] : basename(rstrip(output_dir, '/'))
    calibration_report_pdf(output_dir; label)
end
