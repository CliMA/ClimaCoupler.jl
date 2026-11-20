# Check that a warmup build produced usable images.
#
#     julia --project=experiments/warmup experiments/warmup/check_warmup.jl
#
# Catches three failures that `Pkg.precompile()` does not report:
#
#   * A workload that threw. Each workload catches its own error and records it in
#     `workload_report.txt`, so Pkg sees a successful build.
#   * A link that failed. An image over ~2 GiB fails to link after the workload returns,
#     leaving a stale image alongside `status = ok`.
#   * An image that cached almost nothing, detected by the size floors in `MIN_MIB`.
#
# Exits non-zero if any package is unusable, so it can gate a CI step. Measurements behind
# these checks are in `precompile_info.md`.

const PACKAGES = [
    "AtmosWarmupSetup",
    "AtmosWarmupExplicit",
    "AtmosWarmupStep",
    "OceanWarmup",
    "LandWarmup",
    "SeaIceWarmup",
    "CMIPWarmup",
]

# Minimum image size per package: low enough to survive ordinary version drift, high enough
# to catch an image that cached nothing. Packages high in the DAG cache only their increment
# over the packages beneath them, so their floors are small.
const MIN_MIB = Dict(
    "AtmosWarmupSetup" => 300,
    "AtmosWarmupExplicit" => 150,
    "AtmosWarmupStep" => 150,
    "OceanWarmup" => 200,
    "LandWarmup" => 200,
    "SeaIceWarmup" => 25,
    "CMIPWarmup" => 75,
)

"""
    image_path(pkg)

Return the newest package image for `pkg` in the first depot, or `nothing`.
"""
function image_path(pkg)
    dir = joinpath(DEPOT_PATH[1], "compiled", "v$(VERSION.major).$(VERSION.minor)", pkg)
    isdir(dir) || return nothing
    sos = filter(f -> endswith(f, ".so"), readdir(dir; join = true))
    isempty(sos) ? nothing : argmax(mtime, sos)
end

"""
    source_mtime(pkg)

Return the latest modification time across the sources of `pkg` and the warmup configs.
"""
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
