using Documenter
using CouplerMachine

makedocs(
    sitename = "CouplerMachine",
    format = Documenter.HTML(),
    modules = [CouplerMachine]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "<github.com/CliMA/CouplerMachine.jl.git>",
    push_preview = true,
)
