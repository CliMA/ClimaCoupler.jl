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

#   * An image that built cleanly and cached almost nothing. A workload can run, report
#     success and still produce a near-empty image if the mechanism that serializes
#     *externally-owned* code is not active -- which is what PrecompileTools 1.2 does on
#     Julia >= 1.12. Observed: AtmosWarmupExplicit at 2 MiB where it should be ~1200 MiB,
#     with `status = ok`. Hence the minimum sizes below.

const PACKAGES = [
    "AtmosWarmupSetup",
    "AtmosWarmupExplicit",
    "AtmosWarmupStep",
    "OceanWarmup",
    "LandWarmup",
    "SeaIceWarmup",
    "CMIPWarmup",
]

# A floor, not an expectation: enough to catch an image that cached nothing, loose enough
# not to need updating for ordinary drift. Sizes move legitimately with Julia and
# dependency versions -- the whole stack shrank about 2.7x between Julia 1.11/CUDA 5.11 and
# Julia 1.13/CUDA 6.4, even though the latter also persists GPU kernels into the image
# (measured at ~89 KiB each).
#
# A package that sits high in the DAG is *expected* to be small: it caches only its own
# increment over everything beneath it. SeaIceWarmup is 80 MiB because OceanWarmup carries
# the ocean, and CMIPWarmup is 219 MiB because the six packages below it carry the
# components -- a run against that stack still reached its first step with 5.7 s of setup
# compilation, so the small image reflects good layering, not a failed workload.
const MIN_MIB = Dict(
    "AtmosWarmupSetup" => 300,
    "AtmosWarmupExplicit" => 150,
    "AtmosWarmupStep" => 150,
    "OceanWarmup" => 200,
    "LandWarmup" => 200,
    "SeaIceWarmup" => 25,
    "CMIPWarmup" => 75,
)

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
    floor_mib = get(MIN_MIB, pkg, 0)
    big_enough = mib >= floor_mib
    println(
        rpad(pkg, 24),
        rpad("$(Int(mib)) MiB", 12),
        rpad(fresh ? "yes" : "STALE", 8),
        status,
    )
    fresh ||
        push!(failures, "$pkg: image older than its sources/configs -- link likely failed")
    big_enough || push!(
        failures,
        "$pkg: image is $(Int(mib)) MiB, under the $(floor_mib) MiB floor -- the workload " *
        "ran but cached little or nothing (check PrecompileTools >= 1.3)",
    )
    occursin("status = ok", status) ||
        push!(failures, "$pkg: workload did not succeed -- $status")
end

println()
if isempty(failures)
    println("All warmup packages usable.")
else
    println("UNUSABLE ($(length(failures))):")
    foreach(f -> println("  - ", f), failures)
    println(
        "\nA failure low in the DAG makes the packages above it re-cache its content, which",
    )
    println("can push them over the 2 GiB limit. Fix the lowest failure first and rebuild.")
    exit(1)
end
