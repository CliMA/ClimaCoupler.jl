#=
Add the packages listed in the `UPSTREAM_PACKAGES` environment variable to the
active environment. Entries are separated by spaces and take one of these forms:
- `Name`: the main branch
- `Name@version`: a registered version
- `Name#rev`: a branch, tag, or commit
- `url#rev`: a branch, tag, or commit of the repository at `url`, e.g. a fork

With `--main`, every `Name#rev` and `url#rev` entry uses the main branch of the
registered package instead, for the time-to-first-step baseline.

```bash
UPSTREAM_PACKAGES="ClimaAtmos#abc1234 ClimaCore SurfaceFluxes@0.12" \
    julia --project=experiments/AMIP .buildkite/add_upstream_packages.jl
```
=#

import Pkg

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

use_main = "--main" in ARGS
specs = map(e -> package_spec(e, use_main), split(get(ENV, "UPSTREAM_PACKAGES", "")))
isempty(specs) || Pkg.add(specs)
