# Contributing

Thank you for contributing to `ClimaCoupler`! We encourage Pull Requests (PRs).
Please do not hesitate to ask questions, or to open issues if something seems amiss
or you'd like a new feature.

## Shared developer guidelines

ClimaCoupler.jl vendors the [CliMA DeveloperGuides](https://github.com/CliMA/DeveloperGuides)
as a Git subtree at `docs/dev-guides/`. These shared guides cover cross-repo
standards that apply to all CliMA packages:

- **Architecture & design patterns** — [architectural boundaries](https://github.com/CliMA/DeveloperGuides/blob/main/architecture/architectural_boundaries.md), [cross-repo contracts](https://github.com/CliMA/DeveloperGuides/blob/main/architecture/cross_repo_contracts.md)
- **Code quality** — [code style](https://github.com/CliMA/DeveloperGuides/blob/main/code-quality/code_style.md), [changelogs & versions](https://github.com/CliMA/DeveloperGuides/blob/main/code-quality/changelogs_and_versions.md), [documentation policy](https://github.com/CliMA/DeveloperGuides/blob/main/code-quality/documentation_policy.md)
- **Performance** — GPU patterns, type stability, AD compatibility
- **Workflow** — [PR review](https://github.com/CliMA/DeveloperGuides/blob/main/workflow/review.md), [debugging](https://github.com/CliMA/DeveloperGuides/blob/main/workflow/debugging.md), [CI triage](https://github.com/CliMA/DeveloperGuides/blob/main/workflow/ci_triage.md)

For conventions specific to this repo (directory layout, test groups, key
abstractions), see [`docs/clima_coupler_specific.md`](https://github.com/CliMA/ClimaCoupler.jl/blob/main/docs/clima_coupler_specific.md).

!!! note
    Edits to the shared guides belong in
    [CliMA/DeveloperGuides](https://github.com/CliMA/DeveloperGuides), not in
    the vendored copy here. The subtree is synced monthly via
    `.github/workflows/update_dev_guides.yml`.

## Modern Julia Workflows (general development advice)
For tips on writing, sharing, and optimizing Julia packages — including development workflows,
testing, documentation, and performance — see [Modern Julia Workflows](https://modernjuliaworkflows.org).
The site includes general tips on everything from setting up your editor and managing environments,
to registering your package and reducing compilation latency.

## Other useful tips
- When developing code it's best to work on a branch off of the most recent main.
This can be done by running the following commands, where "initials" corresponds to the first and last initial of the person starting the branch.
```
git checkout main
git pull
git checkout -b initials/branch_name
```

- Make sure you add tests for your code in `test/`, appropriate documentation in `docs/`,
  and descriptive inline comments throughout the code.
  All exported functions and structs must have docstrings.
- When your PR is ready for review, clean up your commit history by squashing to 1 commit per PR
  and make sure your code is current with the main branch by rebasing.

## Continuous integration

After rebasing your branch, you can ask for review. Fill out the template and
provide a clear summary of what your PR does. When a PR is created or
updated, a set of automated tests are run on the PR in our continuous
integration (CI) system.

ClimaCoupler.jl's continuous integration contains both unit tests and integration
tests (coupled simulations).

### Formatting check

The `JuliaFormatter` test checks if the PR is correctly formatted according to
the project's style guidelines. The previous `.dev/climaformat.jl` script has
been discontinued in favor of using the JuliaFormatter package directly.

To format your code, first add JuliaFormatter to your base environment:

```sh
julia -e 'using Pkg; Pkg.add(PackageSpec("JuliaFormatter", v"2.10.1"))'
```

Then, in a Julia REPL, run:

```julia
using JuliaFormatter; format(".")
```


### Documentation

The `Documentation` test rebuilds the documentation for the PR and checks if the docs
are consistent and generate valid output.

To add internal references, for example to another documentation page or API, see the relevant
`Documenter.jl` `@ref` [documentation page](@extref Documenter Named-@refs), example syntax:
`[see contributor guide](@ref "Contributing")` for a page or
`[CoupledSimulation object and constructors](@ref ClimaCoupler.Interfacer.CoupledSimulation)` for an API.

To add external references, for example to another package documentation page or API, see
the documentation of [DocumenterInterLinks](https://juliadocs.org/DocumenterInterLinks.jl/stable/).
Example syntax: `[how to us @ref](@extref Documenter Named-@refs)`.

To add a reference from the literature, see the documentation of
[DocumenterCitations.jl](https://juliadocs.org/DocumenterCitations.jl/stable/), in short, add your reference
to `docs/src/refs.bib` and then refer to it with this syntax `[authors1999](@citet)`.
