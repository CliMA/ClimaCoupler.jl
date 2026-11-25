# Check that the warmup environment resolves to the same package versions as
# experiments/CMIP.
#
#     julia --project=experiments/warmup experiments/warmup/check_version_drift.jl
#
# Exits non-zero if a package in both manifests differs and is not in `ALLOWED`, or if
# either environment is below a version in `FLOORS`.

using TOML

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const MANIFEST = "Manifest-v$(VERSION.major).$(VERSION.minor).toml"

# Packages allowed to differ, with the reason. Only build-time tooling belongs here.
const ALLOWED =
    Dict("PrecompileTools" => "build-time only, absent from the cached code; see FLOORS")

# 1.2 caches no externally-owned code on Julia >= 1.12; see precompile_info.md.
const FLOORS = Dict("PrecompileTools" => (VERSION >= v"1.12" ? v"1.3" : v"0"))

"""
    versions(path)

Return a `Dict` of package name => version string for the manifest at `path`.
"""
function versions(path)
    deps = get(TOML.parsefile(path), "deps", Dict{String, Any}())
    out = Dict{String, String}()
    for (name, entries) in deps, entry in entries
        haskey(entry, "version") && (out[name] = entry["version"])
    end
    out
end

const RUN_MANIFEST = joinpath(REPO, "experiments", "CMIP", MANIFEST)
const WARM_MANIFEST = joinpath(REPO, "experiments", "warmup", MANIFEST)

for p in (RUN_MANIFEST, WARM_MANIFEST)
    isfile(p) || (println("missing manifest: ", p); exit(1))
end

run_versions, warm_versions = versions(RUN_MANIFEST), versions(WARM_MANIFEST)
shared = sort!(collect(intersect(keys(run_versions), keys(warm_versions))))

differing = [
    (n, run_versions[n], warm_versions[n]) for
    n in shared if run_versions[n] != warm_versions[n]
]
mismatched = filter(t -> !haskey(ALLOWED, first(t)), differing)

println("compared $(length(shared)) packages present in both environments")
for (name, runv, warmv) in differing
    haskey(ALLOWED, name) || continue
    println("  allowed: ", name, " CMIP ", runv, ", warmup ", warmv)
    println("    ", ALLOWED[name])
end

# An allowed difference does not exempt either environment from the floor.
for (name, floor) in FLOORS,
    (env, vs) in (("experiments/CMIP", run_versions), ("warmup", warm_versions))

    haskey(vs, name) || continue
    if VersionNumber(vs[name]) < floor
        push!(mismatched, (name, "$env has $(vs[name]), below the $floor floor", ""))
    end
end

if isempty(mismatched)
    println("Warmup and experiments/CMIP agree on every shared package version.")
else
    println("$(length(mismatched)) version mismatch(es):")
    for (name, runv, warmv) in mismatched
        println("  ", rpad(name, 28), "CMIP ", rpad(runv, 12), "warmup ", warmv)
    end
    println()
    println("The warmup images will not be used by a run in experiments/CMIP. Re-seed:")
    println("  cp experiments/CMIP/$MANIFEST experiments/warmup/$MANIFEST")
    println("  julia --project=experiments/warmup -e 'using Pkg; Pkg.resolve()'")
    exit(1)
end
