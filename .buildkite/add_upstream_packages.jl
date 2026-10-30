#=
Add the packages listed in the `UPSTREAM_PACKAGES` environment variable to the
active environment. Entries are separated by spaces and take one of these forms:
- `Name`: the main branch
- `Name@version`: a registered version
- `Name#rev`: a branch, tag, or commit
- `url#rev`: a branch, tag, or commit of the repository at `url`, e.g. a fork

With `--main`, every `Name#rev` and `url#rev` entry uses the main branch of the
registered package instead, for the time-to-first-step baseline.

The `[compat]` entries of these packages are relaxed to `>= 0` in the active
project and in ClimaCoupler's Project.toml.

```bash
UPSTREAM_PACKAGES="ClimaAtmos#abc1234 ClimaCore SurfaceFluxes@0.12" \
    julia --project=experiments/AMIP .buildkite/add_upstream_packages.jl
```
=#

import Pkg
import TOML

package_name(repo) = String(chopsuffix(chopsuffix(basename(repo), ".git"), ".jl"))

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

# Edits the TOML directly because `Pkg.compat` rejects weak dependencies
function relax_compat(project_file, names)
    project = TOML.parsefile(project_file)
    compat = get(project, "compat", Dict{String, Any}())
    relaxed = filter(name -> get(compat, name, ">= 0") != ">= 0", names)
    isempty(relaxed) && return nothing
    @info "Relaxing [compat] entries to >= 0" project_file relaxed = join(relaxed, ", ")
    foreach(name -> compat[name] = ">= 0", relaxed)
    open(io -> TOML.print(io, project; sorted = true), project_file, "w")
    return nothing
end

use_main = "--main" in ARGS
specs = map(e -> package_spec(e, use_main), split(get(ENV, "UPSTREAM_PACKAGES", "")))
names = [@something(spec.name, package_name(spec.url)) for spec in specs]
coupler_project = joinpath(dirname(@__DIR__), "Project.toml")
foreach(p -> relax_compat(p, names), unique([Base.active_project(), coupler_project]))
isempty(specs) || Pkg.add(specs)
