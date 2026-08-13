# Frozen run packages for calibrations.
#
# A package is an immutable snapshot of the repository inside a run directory.
# The calibration executes from the snapshot, so later edits or rebases in the
# working repository cannot change a running or resumed calibration.
#
# Layout of a run root:
#   <run_root>/package/ClimaCoupler.jl/   snapshot: git archive + dirty patch
#   <run_root>/package/provenance.toml    what was launched, from where
#   <run_root>/package/dirty.patch        uncommitted tracked changes at launch
#   <run_root>/package/driver.pbs         the launch and resume entry point
#   <run_root>/package/fingerprint.sha256 per-file hashes of the package
#   <run_root>/output/                    ClimaCalibrate output_dir
#   <run_root>/logs/                      driver log
#
# The freeze works because experiments/AMIP resolves ClimaCoupler through the
# relative path "../.." in its Manifest, and every other dependency is pinned
# by content hash to an immutable directory in the depot. Activating the
# snapshot's experiments/AMIP therefore runs the snapshot's code.

import Dates
import SHA
import TOML

const PACKAGING_JULIA_BIN = "/glade/campaign/univ/ucit0011/software/julia/julia-1.11.3/bin/julia"

sha256_file(path) = bytes2hex(SHA.sha256(read(path)))

snapshot_root(run_root) = joinpath(run_root, "package", "ClimaCoupler.jl")
snapshot_env(run_root) = joinpath(snapshot_root(run_root), "experiments", "AMIP")
provenance_path(run_root) = joinpath(run_root, "package", "provenance.toml")
fingerprint_path(run_root) = joinpath(run_root, "package", "fingerprint.sha256")

"""
    create_snapshot(repo_root, run_root)

Copy the repository into `<run_root>/package/ClimaCoupler.jl` with
`git archive` at HEAD, then apply the uncommitted tracked changes as a patch.
Untracked files are never captured. Returns a NamedTuple with the git state.
"""
function create_snapshot(repo_root, run_root)
    dest = snapshot_root(run_root)
    mkpath(dest)

    sha = readchomp(`git -C $repo_root rev-parse HEAD`)
    branch = readchomp(`git -C $repo_root rev-parse --abbrev-ref HEAD`)
    describe = readchomp(`git -C $repo_root describe --always --tags`)

    patch_path = joinpath(run_root, "package", "dirty.patch")
    open(patch_path, "w") do io
        run(pipeline(`git -C $repo_root diff HEAD`; stdout = io))
    end
    dirty = filesize(patch_path) > 0

    status = readchomp(`git -C $repo_root status --porcelain`)
    untracked = [
        strip(line[4:end]) for
        line in split(status, '\n'; keepempty = false) if startswith(line, "??")
    ]

    run(pipeline(`git -C $repo_root archive $sha`, `tar -x -C $dest`))
    dirty && run(Cmd(`git apply $patch_path`; dir = dest))

    return (;
        sha,
        branch,
        describe,
        dirty,
        dirty_patch_sha256 = dirty ? sha256_file(patch_path) : "",
        untracked,
    )
end

"""
    output_dir_name_from_config(config_path)

Read the repository-relative output directory name from a config file. Config
files set `output_dir = joinpath(pkgdir(ClimaCoupler), "<name>")`, and under
the snapshot environment that path lands inside the snapshot. The launcher
places a symlink there that points at `<run_root>/output`.
"""
function output_dir_name_from_config(config_path)
    text = read(config_path, String)
    m = match(r"output_dir\s*=\s*joinpath\(pkgdir\(ClimaCoupler\),\s*\"([^\"]+)\"\)", text)
    isnothing(m) && error(
        "Cannot find `output_dir = joinpath(pkgdir(ClimaCoupler), \"...\")` in $config_path",
    )
    return String(m.captures[1])
end

"""
    link_output!(run_root, config_path)

Create `<run_root>/output` and symlink the snapshot's computed output
directory to it.
"""
function link_output!(run_root, config_path)
    outdir = joinpath(run_root, "output")
    mkpath(outdir)
    name = output_dir_name_from_config(config_path)
    linkpath = joinpath(snapshot_root(run_root), name)
    islink(linkpath) || symlink(outdir, linkpath)
    return outdir
end

"""
    validate_snapshot(run_root, config_rel)

Assert that every file the run needs exists inside the snapshot: the config,
the coupler YAML it names, the atmos YAML, and each coupler TOML. This catches
untracked config files at launch instead of hours into a job.
"""
function validate_snapshot(run_root, config_rel)
    snap = snapshot_root(run_root)
    config_path = joinpath(snap, config_rel)
    isfile(config_path) || error(
        "Config $config_rel is not in the snapshot. Untracked files are not " *
        "captured. Commit the config first.",
    )
    text = read(config_path, String)
    m = match(r"\"config\",\s*\"amip_configs\",\s*\"([^\"]+)\"", text)
    isnothing(m) && return nothing
    yaml_path = joinpath(snap, "config", "amip_configs", String(m.captures[1]))
    isfile(yaml_path) || error("Coupler YAML $(m.captures[1]) is not in the snapshot")
    yaml = read(yaml_path, String)
    for pat in (r"atmos_config_file:\s*\"([^\"]+)\"", r"coupler_toml:\s*\[\"([^\"]+)\"\])?")
        for mm in eachmatch(pat, yaml)
            rel = String(mm.captures[1])
            isfile(joinpath(snap, rel)) ||
                error("Snapshot is missing $rel referenced by $(m.captures[1])")
        end
    end
    return nothing
end

"""
    instantiate_snapshot(run_root; julia_bin = PACKAGING_JULIA_BIN)

Run `Pkg.instantiate()` in the snapshot environment. The Manifest is already
resolved, so this verifies that the depot can serve every pinned dependency
while the live and snapshot environments are still identical.
"""
function instantiate_snapshot(run_root; julia_bin = PACKAGING_JULIA_BIN)
    env = snapshot_env(run_root)
    run(
        Cmd(
            `$julia_bin --project=$env -e 'using Pkg; Pkg.instantiate()'`;
            dir = snapshot_root(run_root),
        ),
    )
    return nothing
end

"Key dependency versions from the snapshot Manifest, for the provenance file."
function key_deps(run_root)
    manifest = joinpath(snapshot_env(run_root), "Manifest-v1.11.toml")
    parsed = TOML.parsefile(manifest)
    deps = get(parsed, "deps", Dict())
    out = Dict{String, Any}()
    for name in ("ClimaAtmos", "ClimaCore", "ClimaCalibrate", "ClimaLand", "RRTMGP")
        haskey(deps, name) || continue
        entry = deps[name][1]
        out[name] = Dict(
            k => entry[k] for
            k in ("version", "git-tree-sha1", "repo-rev") if haskey(entry, k)
        )
    end
    return out
end

"""
    write_provenance(run_root; git_state, config_rel, extra = Dict())

Write `<run_root>/package/provenance.toml`.
"""
function write_provenance(run_root; git_state, config_rel, extra = Dict())
    snap = snapshot_root(run_root)
    manifest = joinpath(snapshot_env(run_root), "Manifest-v1.11.toml")
    doc = Dict(
        "schema_version" => 1,
        "launch" => Dict(
            "time" => string(Dates.now()),
            "user" => get(ENV, "USER", "unknown"),
            "host" => gethostname(),
            "run_root" => run_root,
            "config" => config_rel,
        ),
        "git" => Dict(
            "sha" => git_state.sha,
            "branch" => git_state.branch,
            "describe" => git_state.describe,
            "dirty" => git_state.dirty,
            "dirty_patch_sha256" => git_state.dirty_patch_sha256,
            "untracked_ignored" => git_state.untracked,
        ),
        "julia" => Dict("version" => string(VERSION), "binary" => PACKAGING_JULIA_BIN),
        "environment" => Dict(
            "manifest_sha256" => sha256_file(manifest),
            "project_sha256" =>
                sha256_file(joinpath(snapshot_env(run_root), "Project.toml")),
            "key_deps" => key_deps(run_root),
        ),
        "config_sha256" => sha256_file(joinpath(snap, config_rel)),
    )
    merge!(doc, extra)
    open(provenance_path(run_root), "w") do io
        TOML.print(io, doc)
    end
    return nothing
end

"""
    fingerprint_package(run_root)

Hash every regular file under `<run_root>/package` except the fingerprint file
itself and write the sorted `hash  relative/path` lines to
`package/fingerprint.sha256`.
"""
function fingerprint_package(run_root)
    pkg = joinpath(run_root, "package")
    fp = fingerprint_path(run_root)
    lines = String[]
    for (root, _, files) in walkdir(pkg; follow_symlinks = false)
        for f in files
            path = joinpath(root, f)
            islink(path) && continue
            path == fp && continue
            rel = relpath(path, pkg)
            push!(lines, sha256_file(path) * "  " * rel)
        end
    end
    sort!(lines)
    write(fp, join(lines, "\n") * "\n")
    return nothing
end

"""
    verify_package(run_root; strict = true)

Recompute the package fingerprint and compare it to the recorded one. Changed
or missing files are an error when `strict` is true, because a resumed run must
execute exactly the code that was launched. Files that appeared after launch
(for example worker logs written into the snapshot working directory) only
produce a warning. Also warns when a Manifest-pinned dependency is missing from
the depot, which happens after `Pkg.gc()`.
"""
function verify_package(run_root; strict = true)
    pkg = joinpath(run_root, "package")
    fp = fingerprint_path(run_root)
    isfile(fp) || error("No fingerprint at $fp. Not a packaged run.")
    recorded = Dict{String, String}()
    for line in eachline(fp)
        h, rel = split(line, "  "; limit = 2)
        recorded[String(rel)] = String(h)
    end

    changed = String[]
    missing_files = String[]
    for (rel, h) in recorded
        path = joinpath(pkg, rel)
        if !isfile(path)
            push!(missing_files, rel)
        elseif sha256_file(path) != h
            push!(changed, rel)
        end
    end
    extra = String[]
    for (root, _, files) in walkdir(pkg; follow_symlinks = false)
        for f in files
            path = joinpath(root, f)
            islink(path) && continue
            path == fp && continue
            rel = relpath(path, pkg)
            haskey(recorded, rel) || push!(extra, rel)
        end
    end

    isempty(extra) ||
        @warn "Files appeared in the package after launch" n = length(extra) sample =
            first(extra, 5)
    if !(isempty(changed) && isempty(missing_files))
        msg = "Package fingerprint mismatch. Changed: $(changed). Missing: $(missing_files)."
        strict ? error(msg) : @warn msg
    end

    # Depot check: every pinned dependency must still exist.
    manifest = joinpath(snapshot_env(run_root), "Manifest-v1.11.toml")
    if isfile(manifest)
        parsed = TOML.parsefile(manifest)
        depot = joinpath(first(DEPOT_PATH), "packages")
        gone = String[]
        for (name, entries) in get(parsed, "deps", Dict())
            entry = entries[1]
            haskey(entry, "git-tree-sha1") || continue
            slugs =
                readdir(joinpath(depot, name); join = false, sort = false) |>
                x -> isdir(joinpath(depot, name)) ? x : String[]
            isempty(slugs) && push!(gone, name)
        end
        isempty(gone) ||
            @warn "Depot no longer holds these pinned packages. Pkg.gc may have removed them" gone
    end
    return true
end

"""
    find_partial_members(output_dir)

Return member directories whose checkpoint says started but not completed. A
member resumed from a mid-run model checkpoint writes wrong-dated monthly
diagnostics and crashes the observation map, so a resume must clean these
first.
"""
function find_partial_members(output_dir)
    partial = String[]
    isdir(output_dir) || return partial
    for it in filter(d -> startswith(d, "iteration_"), readdir(output_dir))
        itdir = joinpath(output_dir, it)
        isdir(itdir) || continue
        for mem in filter(d -> startswith(d, "member_"), readdir(itdir))
            ckpt = joinpath(itdir, mem, "checkpoint.txt")
            isfile(ckpt) || continue
            strip(read(ckpt, String)) == "completed" && continue
            push!(partial, joinpath(itdir, mem))
        end
    end
    return partial
end

"""
    clean_partial_members!(output_dir)

Reset partially-run members so they rerun from the start of the iteration.
Keeps `parameters.toml`, removes the checkpoint and the model output.
"""
function clean_partial_members!(output_dir)
    for mdir in find_partial_members(output_dir)
        for entry in readdir(mdir)
            entry == "parameters.toml" && continue
            rm(joinpath(mdir, entry); recursive = true, force = true)
        end
        @info "Reset partial member" mdir
    end
    return nothing
end
