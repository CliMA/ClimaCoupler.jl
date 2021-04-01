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
#=deploydocs(
    repo = "<https://github.com/CliMA/CouplerMachine>",
    push_preview = true,
)=#
