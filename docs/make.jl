using Documenter
using CouplerMachine
using Literate

const EXPERIMENTS_DIR = joinpath(@__DIR__, "..", "experiments")
const OUTPUT_DIR      = joinpath(@__DIR__, "..", "docs/src/generated")

experiments = [
               "DesignTests/simple_2testcomp.jl",
              ]

for experiment in experiments
    experiment_filepath = joinpath(EXPERIMENTS_DIR, experiment)
    Literate.markdown(experiment_filepath, OUTPUT_DIR, documenter=true)
end

experiment_pages = [
                    "Simple Two Component Test" => "generated/simple_2testcomp.md",
                   ]

pages = Any[
    "Home" => "index.md",
    "Examples" => experiment_pages,
    "Coupler Interface" => ["couplerstate.md", "timestepping.md"],
]


makedocs(
    sitename = "CouplerMachine",
    format = Documenter.HTML(),
    modules = [CouplerMachine],
    pages = pages,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "<github.com/CliMA/CouplerMachine.jl.git>",
    devbranch = "main",
    push_preview = true,
)
