#=
Add the packages listed in the `UPSTREAM_PACKAGES` environment variable to the
active environment. Entries are separated by spaces and take one of these forms:
- `Name`: the main branch
- `Name@version`: a registered version
- `Name#rev`: a branch, tag, or commit
- `url#rev`: a branch, tag, or commit of the repository at `url`, e.g. a fork

With `--main`, every `Name#rev` and `url#rev` entry uses the main branch of the
registered package instead, for the time-to-first-step baseline.

With `--skip-missing`, packages that are not in the active environment are
skipped.

With `--relax-compat`, the `[compat]` entries of the added packages in the
active project and in ClimaCoupler become `>= 0`, so versions outside those
bounds resolve.

```bash
UPSTREAM_PACKAGES="ClimaAtmos#abc1234 ClimaCore SurfaceFluxes@0.12" \
    julia --project=experiments/AMIP .buildkite/add_upstream_packages.jl
```
=#

import Pkg
import TOML

package_name(repo) = String(chopsuffix(chopsuffix(basename(repo), ".git"), ".jl"))
spec_name(spec) = isnothing(spec.name) ? package_name(spec.repo.source) : spec.name

function package_spec(entry, use_main)
    if occursin('@', entry) && !occursin('/', entry)
        name, version = split(entry, '@'; limit = 2)
        return Pkg.PackageSpec(; name = package_name(name), version = String(version))
    end
    repo, rev = occursin('#', entry) ? split(entry, '#'; limit = 2) : (entry, "main")
    use_main && return Pkg.PackageSpec(; name = package_name(repo), rev = "main")
    occursin('/', repo) && return Pkg.PackageSpec(; url = String(repo), rev = String(rev))
    return Pkg.PackageSpec(; name = package_name(repo), rev = String(rev))
end

# Edit the TOML directly since `Pkg.compat` rejects weakdeps. Keep the entries (as
# `>= 0`) since `Aqua.test_deps_compat` requires them.
function relax_compat!(project_file, names)
    project = TOML.parsefile(project_file)
    compat = get(project, "compat", Dict{String, Any}())
    relaxed = filter(name -> haskey(compat, name), names)
    isempty(relaxed) && return nothing
    @info "Relaxing [compat] entries" project_file relaxed
    for name in relaxed
        compat[name] = ">= 0"
    end
    open(io -> TOML.print(io, project), project_file, "w")
    return nothing
end

use_main = "--main" in ARGS
specs = map(e -> package_spec(e, use_main), split(get(ENV, "UPSTREAM_PACKAGES", "")))

if "--skip-missing" in ARGS
    env_names = Set(dep.name for dep in values(Pkg.dependencies()))
    for spec in filter(spec -> spec_name(spec) ∉ env_names, specs)
        @info "Skipping $(spec_name(spec)), which is not in the environment"
    end
    filter!(spec -> spec_name(spec) in env_names, specs)
end

if "--relax-compat" in ARGS
    climacoupler_project = joinpath(dirname(@__DIR__), "Project.toml")
    for project_file in (Base.active_project(), climacoupler_project)
        relax_compat!(project_file, map(spec_name, specs))
    end
end

isempty(specs) || Pkg.add(specs)
