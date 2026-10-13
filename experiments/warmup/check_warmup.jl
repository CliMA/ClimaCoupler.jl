# Verify that a warmup build actually produced usable images.
#
#     julia --project=experiments/warmup experiments/warmup/check_warmup.jl
#
# Two failure modes are invisible to the workloads themselves and need checking here:
#
#   * A workload that threw. Each package catches its own exception so that one broken
#     workload cannot break precompilation of the whole environment, and records the result
#     in `workload_report.txt`. Pkg therefore sees a successful precompile and happily
#     builds everything above it -- so a failure low in the DAG silently degrades the
#     packages that depend on it. (Set `WARMUP_STRICT=1` to rethrow instead, which makes Pkg
#     skip dependents.)
#
#   * A link that failed. Package images cannot exceed ~2 GiB: beyond that, 32-bit ELF
#     relocations cannot encode the distance into the serialized blob and `lld` refuses.
#     This happens *after* the workload returns, so no workload-level error handling sees
#     it. The symptom is a stale image alongside `status = ok`.
#
# Exits non-zero if any package is unusable, so it can gate a CI step.

const PACKAGES = [
    "AtmosWarmupSetup",
    "AtmosWarmupExplicit",
    "AtmosWarmupStep",
    "OceanWarmup",
    "LandWarmup",
    "SeaIceWarmup",
    "CMIPWarmup",
]

"Newest package image for `pkg`, or `nothing`."
function image_path(pkg)
    dir = joinpath(DEPOT_PATH[1], "compiled", "v$(VERSION.major).$(VERSION.minor)", pkg)
    isdir(dir) || return nothing
    sos = filter(f -> endswith(f, ".so"), readdir(dir; join = true))
    isempty(sos) ? nothing : argmax(mtime, sos)
end

"Most recent modification time across a package's sources and the configs it reads."
function source_mtime(pkg)
    here = @__DIR__
    srcs = String[]
    pkgdir = joinpath(here, pkg)
    isdir(pkgdir) && for (root, _, files) in walkdir(pkgdir)
        append!(srcs, joinpath(root, f) for f in files if endswith(f, ".jl"))
    end
    cfgdir = joinpath(here, "configs")
    isdir(cfgdir) && append!(srcs, readdir(cfgdir; join = true))
    isempty(srcs) ? 0.0 : maximum(mtime, srcs)
end

failures = String[]
println(rpad("package", 24), rpad("image", 12), rpad("fresh", 8), "workload")
println("-"^72)
for pkg in PACKAGES
    so = image_path(pkg)
    report = joinpath(@__DIR__, pkg, "workload_report.txt")
    status = isfile(report) ? strip(first(eachline(report))) : "(no report)"
    if so === nothing
        println(rpad(pkg, 24), rpad("MISSING", 12), rpad("-", 8), status)
        push!(failures, "$pkg: no image (did the link fail? check for `out of range`)")
        continue
    end
    mib = round(filesize(so) / 2^20, digits = 0)
    # A stale image means the link failed after the workload ran: the old image survives.
    fresh = mtime(so) >= source_mtime(pkg)
    println(rpad(pkg, 24), rpad("$(Int(mib)) MiB", 12), rpad(fresh ? "yes" : "STALE", 8), status)
    fresh || push!(failures, "$pkg: image older than its sources/configs -- link likely failed")
    occursin("status = ok", status) ||
        push!(failures, "$pkg: workload did not succeed -- $status")
end

println()
if isempty(failures)
    println("All warmup packages usable.")
else
    println("UNUSABLE ($(length(failures))):")
    foreach(f -> println("  - ", f), failures)
    println("\nA failure low in the DAG makes the packages above it re-cache its content, which")
    println("can push them over the 2 GiB limit. Fix the lowest failure first and rebuild.")
    exit(1)
end
