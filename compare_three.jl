# Three-way, time-resolved comparison.
#
#   A = cc_serialA, B = cc_serialB  -> noise floor (identical configs)
#   A vs C = cc_concurrent          -> the signal
#
# Comparing whole files is misleading here: the runs diverge chaotically from
# roundoff, so late times have a large floor while early times do not. We
# therefore compare slice by slice in time and ask, at each time, whether the
# concurrent difference is comparable to the serial-vs-serial difference.

import NCDatasets, JLD2
import Printf: @printf

const RUN_A = length(ARGS) > 0 ? ARGS[1] : "cc_serialA"
const RUN_B = length(ARGS) > 1 ? ARGS[2] : "cc_serialB"
const RUN_C = length(ARGS) > 2 ? ARGS[3] : "cc_concurrent"
# accept "name" (uses output_active) or "name/output_0000" (explicit run dir)
resolve_run(n) =
    occursin('/', n) ? joinpath("output", n) : joinpath("output", n, "output_active")
const A = resolve_run(RUN_A)
const B = resolve_run(RUN_B)
const C = resolve_run(RUN_C)
const COMPONENTS =
    ("clima_atmos", "clima_coupler", "clima_land", "clima_ocean", "clima_seaice")
const COORDS = Set(["lat", "lon", "time", "time_bnds", "z", "z_bnds", "date", "long"])

maxabsdiff(x, y) = begin
    ok = isfinite.(x) .& isfinite.(y)
    any(ok) || return 0.0
    maximum(abs.(Float64.(x[ok]) .- Float64.(y[ok])))
end

"time-sliced max-abs differences: returns vector over time of (ab, ac)"
function nc_slices(fn, comp)
    out = Tuple{String, Any, Float64, Float64}[]   # var, time, ab, ac
    pa, pb, pc = joinpath(A, comp, fn), joinpath(B, comp, fn), joinpath(C, comp, fn)
    (isfile(pa) && isfile(pb) && isfile(pc)) || return out
    NCDatasets.NCDataset(pa) do da
        NCDatasets.NCDataset(pb) do db
            NCDatasets.NCDataset(pc) do dc
                times = haskey(da, "time") ? Array(da["time"]) : [nothing]
                for vn in keys(da)
                    vn in COORDS && continue
                    haskey(db, vn) && haskey(dc, vn) || continue
                    a = Array(da[vn])
                    eltype(a) <: Number || continue
                    b = Array(db[vn]);
                    c = Array(dc[vn])
                    size(a) == size(b) == size(c) || continue
                    # Locate the time axis by NAME -- in these files it is the
                    # FIRST dimension, not the last.
                    dn = collect(NCDatasets.dimnames(da[vn]))
                    tdim = findfirst(==("time"), dn)
                    if tdim !== nothing &&
                       length(times) > 1 &&
                       size(a, tdim) == length(times)
                        for i in 1:length(times)
                            sa = selectdim(a, tdim, i);
                            sb = selectdim(b, tdim, i);
                            sc = selectdim(c, tdim, i)
                            push!(
                                out,
                                (vn, times[i], maxabsdiff(sa, sb), maxabsdiff(sa, sc)),
                            )
                        end
                    else
                        push!(out, (vn, nothing, maxabsdiff(a, b), maxabsdiff(a, c)))
                    end
                end
            end
        end
    end
    return out
end

"Oceananigans JLD2: timeseries/<field>/<iteration>"
function jld2_slices(fn, comp)
    out = Tuple{String, Any, Float64, Float64}[]
    pa, pb, pc = joinpath(A, comp, fn), joinpath(B, comp, fn), joinpath(C, comp, fn)
    (isfile(pa) && isfile(pb) && isfile(pc)) || return out
    JLD2.jldopen(pa) do fa
        JLD2.jldopen(pb) do fb
            JLD2.jldopen(pc) do fc
                haskey(fa, "timeseries") || return
                for field in keys(fa["timeseries"])
                    for iter in keys(fa["timeseries"][field])
                        k = "timeseries/$field/$iter"
                        local a, b, c
                        try
                            a = fa[k];
                            b = fb[k];
                            c = fc[k]
                        catch
                            continue
                        end
                        (a isa AbstractArray && eltype(a) <: Number) || continue
                        size(a) == size(b) == size(c) || continue
                        push!(out, (field, iter, maxabsdiff(a, b), maxabsdiff(a, c)))
                    end
                end
            end
        end
    end
    return out
end

function main()
    allrows = Tuple{String, String, Any, Float64, Float64}[]  # file, var, time, ab, ac
    for comp in COMPONENTS
        d = joinpath(A, comp)
        isdir(d) || continue
        for fn in sort(readdir(d))
            rows =
                endswith(fn, ".nc") ? nc_slices(fn, comp) :
                endswith(fn, ".jld2") ? jld2_slices(fn, comp) : continue
            for (vn, t, ab, ac) in rows
                push!(allrows, ("$comp/$fn", vn, t, ab, ac))
            end
        end
    end

    println("="^86)
    println("THREE-WAY COMPARISON   (A=$RUN_A  B=$RUN_B  C=$RUN_C)")
    println("  A-B = noise floor between two identical serial runs")
    println("  A-C = serial vs concurrent")
    println("="^86)

    # group by time value so the growth of the floor is visible
    bytime = Dict{Any, Vector{Tuple{String, String, Float64, Float64}}}()
    for (f, v, t, ab, ac) in allrows
        push!(get!(bytime, t, []), (f, v, ab, ac))
    end

    @printf(
        "\n%-14s %8s %14s %14s %10s\n",
        "time",
        "fields",
        "worst A-B",
        "worst A-C",
        "ratio"
    )
    # times are Float64 (NetCDF) or String iteration numbers (Oceananigans JLD2)
    sortkey(x) =
        x === nothing ? -Inf :
        x isa AbstractString ? (something(tryparse(Float64, x), Inf)) : Float64(x)
    for t in sort(collect(keys(bytime)), by = sortkey)
        rows = bytime[t]
        wab = maximum(r[3] for r in rows)
        wac = maximum(r[4] for r in rows)
        ratio = wab > 0 ? wac / wab : (wac > 0 ? Inf : 1.0)
        @printf(
            "%-14s %8d %14.5g %14.5g %10.3g\n",
            t === nothing ? "static" : string(t),
            length(rows),
            wab,
            wac,
            ratio
        )
    end

    # any field where concurrent is far outside the serial envelope
    println("\n", "-"^86)
    println("Fields where A-C exceeds A-B by more than 10x (and A-C is not negligible):")
    flagged = [
        (f, v, t, ab, ac) for
        (f, v, t, ab, ac) in allrows if ac > 10 * max(ab, 0) && ac > 1e-10
    ]
    if isempty(flagged)
        println("  none")
    else
        sort!(flagged, by = x -> -x[5])
        for (f, v, t, ab, ac) in first(flagged, 15)
            @printf(
                "  %-46s %-14s t=%-10s A-B=%.4g  A-C=%.4g\n",
                f,
                v,
                t === nothing ? "static" : string(t),
                ab,
                ac
            )
        end
        println("  (", length(flagged), " total)")
    end

    wab_all = maximum(r[4] for r in allrows)
    wac_all = maximum(r[5] for r in allrows)
    println("\n", "="^86)
    @printf("overall worst A-B = %.6g\noverall worst A-C = %.6g\n", wab_all, wac_all)
    println(
        isempty(flagged) ?
        "VERDICT: concurrent stays within the serial-vs-serial envelope everywhere." :
        "VERDICT: $(length(flagged)) field/time slices exceed the envelope -- inspect above.",
    )
    println("="^86)
end

main()
