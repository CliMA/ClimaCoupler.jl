# Verify the warmup environment resolves to the same package versions as experiments/CMIP.
#
#     julia --project=experiments/warmup experiments/warmup/check_version_drift.jl
#
# A cached specialization is keyed on the package versions it was built from. If this
# environment resolves ClimaCore 1.0.8 while a run resolves 1.0.5, the images are simply
# not applicable and the run recompiles -- the same silent failure that
# `check_config_drift.jl` guards for model selections.
#
# This is not hypothetical: a fresh resolve under Julia 1.13 floated to ClimaCore 1.0.8,
# ClimaAtmos 0.42.18 and ClimaLand 1.15.0 against CMIP's 1.0.5 / 0.42.17 / 1.14.1. The
# manifest here is therefore seeded from CMIP's and only the warmup-only extras are
# resolved on top; this script checks that the seeding has not since come apart.
#
# Exits non-zero if any package present in both manifests resolves to a different version.

using TOML

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const MANIFEST = "Manifest-v$(VERSION.major).$(VERSION.minor).toml"

# Packages allowed to differ, with the reason. Only build-time tooling belongs here: a
# science package that differs invalidates the cache, which is the whole point of this
# check. Anything listed is reported but does not fail the run.
const ALLOWED = Dict(
    "PrecompileTools" =>
        "1.2 does not cache externally-owned code on Julia >= 1.12 (measured: 0.02 MiB " *
        "vs 143.55 MiB for the same workload), so this environment needs >= 1.3 even " *
        "though experiments/CMIP still pins 1.2.1. It is a build-time tool and does not " *
        "appear in the cached code.",
)

"Map package name => version string for every versioned entry in a manifest."
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

differing = [(n, run_versions[n], warm_versions[n]) for n in shared if run_versions[n] != warm_versions[n]]
mismatched = filter(t -> !haskey(ALLOWED, first(t)), differing)

println("compared $(length(shared)) packages present in both environments")
for (name, runv, warmv) in differing
    haskey(ALLOWED, name) || continue
    println("  allowed: ", name, " CMIP ", runv, ", warmup ", warmv)
    println("    ", ALLOWED[name])
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
