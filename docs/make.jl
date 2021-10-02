using Documenter
using CouplerMachine
using Literate
using Pkg

const EXPERIMENTS_DIR = joinpath(@__DIR__, "..", "experiments")
const OUTPUT_DIR      = joinpath(@__DIR__, "..", "docs/src/generated")

# experiments = [
#               ]

# for experiment in experiments
#     experiment_filepath = joinpath(EXPERIMENTS_DIR, experiment)
#     Literate.markdown(experiment_filepath, OUTPUT_DIR, documenter=true)
# end

# "Page Name" => "path/to/generated_file.md"

# tutorial
# - generate tutorial files
coupler_dir =  dirname(dirname(pathof(CouplerMachine)));
tutorial_dir = string(coupler_dir,"/experiments/ClimaCore/tc1_heat-diffusion-with-slab/")
tutorial_name = "run"
Pkg.activate(string(tutorial_dir,"/Project.toml"))
Literate.markdown(
                string(tutorial_dir, tutorial_name, ".jl");
                execute = true,
                documenter = false,
                )

# - move tutorial files in src
mv(string(coupler_dir, "/docs/", tutorial_name, ".md"), string(coupler_dir,"/docs/src/",tutorial_name, ".md"), force = true)
files = readdir(string(coupler_dir,"/docs/images/"))
png_files = filter(endswith(".png"),files)
for (i,file) in enumerate(png_files)
    mv(string(coupler_dir, "/docs/images/", file), string(coupler_dir,"/docs/src/images/",file), force = true)
end
rm(string(coupler_dir, "/docs/images/"), force = true)

# pages layout
experiment_pages = [ "run.md" ]
interface_pages = ["couplerstate.md", "timestepping.md"]

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
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
