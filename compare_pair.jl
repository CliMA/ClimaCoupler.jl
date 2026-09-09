# Numerically compare the diagnostic output of two runs.
#
#   julia --project=experiments/CMIP compare_pair.jl cc_serialA cc_serialB
#
# Byte comparison is useless here: NetCDF embeds creation timestamps, so files
# from two identical runs always differ on disk. This compares the numbers.

import NCDatasets, JLD2

include("compare_runs_lib.jl")

function compare_pair(runA, runB; label = "")
    # accept "name" (uses output_active) or "name/output_0000" (explicit run dir)
    resolve(n) = occursin('/', n) ? joinpath("output", n) : joinpath("output", n, "output_active")
    A = resolve(runA)
    B = resolve(runB)
    worst = 0.0
    worstwhere = ""
    nfields = 0
    nfiles = 0
    changed = Tuple{String, Float64}[]

    for comp in ("clima_atmos", "clima_coupler", "clima_land", "clima_ocean", "clima_seaice")
        d = joinpath(A, comp)
        isdir(d) || continue
        for fn in sort(readdir(d))
            (endswith(fn, ".nc") || endswith(fn, ".jld2")) || continue
            pa = joinpath(A, comp, fn)
            pb = joinpath(B, comp, fn)
            isfile(pb) || continue
            da, db = load_any(pa), load_any(pb)
            nfiles += 1
            filemax = 0.0
            for k in intersect(keys(da), keys(db))
                s = diffstats(da[k], db[k])
                s[1] === nothing && continue
                nfields += 1
                s[1] > filemax && (filemax = s[1])
                if s[1] > worst
                    worst = s[1]
                    worstwhere = "$comp/$fn : $k"
                end
            end
            filemax > 0 && push!(changed, ("$comp/$fn", filemax))
        end
    end

    println("-"^72)
    println("$label  $runA  vs  $runB")
    println("  files compared        : $nfiles")
    println("  numeric fields        : $nfields")
    println("  files differing       : $(length(changed)) / $nfiles")
    println("  worst abs difference  : $worst")
    isempty(worstwhere) || println("  worst at              : $worstwhere")
    if !isempty(changed)
        sort!(changed, by = x -> -x[2])
        println("  largest differences:")
        for (f, v) in first(changed, 8)
            println("    ", rpad(f, 46), v)
        end
    end
    return worst
end

if abspath(PROGRAM_FILE) == @__FILE__
    a, b = ARGS[1], ARGS[2]
    label = length(ARGS) > 2 ? ARGS[3] : ""
    compare_pair(a, b; label)
end
