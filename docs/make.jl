using Documenter
using CouplerMachine
using Literate

const EXPERIMENTS_DIR = joinpath(@__DIR__, "..", "experiments")
const OUTPUT_DIR      = joinpath(@__DIR__, "..", "docs/src/generated")

experiments = [
               "ClimateMachine/DesignTests/simple_2testcomp.jl",
               "ClimateMachine/AdvectionDiffusion/run_script_v2.jl"
              ]

for experiment in experiments
    experiment_filepath = joinpath(EXPERIMENTS_DIR, experiment)
    Literate.markdown(experiment_filepath, OUTPUT_DIR, documenter=true)
end

experiment_pages = [
                    "Vertical Column Heat Diffusion" => "generated/simple_2testcomp.md",
                    "Advection-diffusion on a Sphere" => "generated/run_script_v2.md",
                   ]

interface_pages = ["couplerstate.md", "timestepping.md", "coupledmodel.md"]

pages = Any[
    "Home" => "index.md",
    "Examples" => experiment_pages,
    "Coupler Interface" => interface_pages,
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
