using ClimaCoupler
using Documenter, Literate, Pkg

const COUPLER_DIR = joinpath(@__DIR__, "..")
const EXPERIMENTS_DIR = joinpath(@__DIR__, "..", "experiments")
const OUTPUT_DIR = joinpath(@__DIR__, "src/generated")

# tutorials & experiments
# - generate tutorial files:

# sea breeze tutorial
tutorial_dir = joinpath(EXPERIMENTS_DIR, "ClimaCore/sea_breeze/")
tutorial_name = "run.jl"
Pkg.activate(tutorial_dir)
Pkg.instantiate()
include(joinpath(tutorial_dir, "run.jl"))
# Literate.markdown(joinpath(tutorial_dir, tutorial_name), OUTPUT_DIR; execute = true, documenter = false)
# execute Literate on all julia files
tutorial_files = filter(x -> last(x, 3) == ".jl", readdir(tutorial_dir))
map(x -> Literate.markdown(joinpath(tutorial_dir, x), OUTPUT_DIR; execute = true, documenter = false), tutorial_files)

# - move tutorial files to docs/src
# IMAGE_DIR = joinpath(tutorial_dir, "images/")
# files = readdir(IMAGE_DIR)
# png_files = filter(endswith(".png"), files)
# for file in png_files
#     mkpath(joinpath(OUTPUT_DIR, "images/"))
#     cp(joinpath(IMAGE_DIR, file), joinpath(OUTPUT_DIR, "images/", file), force = true)
# end

# pages layout
# experiment_pages = ["Sea Breeze" => ["Overview" => README.md, "Atmosphere" => ]]
# interface_pages = ["couplerstate.md", "timestepping.md"]

# pages = Any["Home" => "index.md", "Examples" => experiment_pages, "Coupler Interface" => interface_pages]


# makedocs(sitename = "ClimaCoupler.jl", format = Documenter.HTML(), modules = [ClimaCoupler], pages = pages)

# # Documenter can also automatically deploy documentation to gh-pages.
# # See "Hosting Documentation" and deploydocs() in the Documenter manual
# # for more information.
# deploydocs(repo = "<github.com/CliMA/ClimaCoupler.jl.git>", push_preview = true, devbranch = "main", forcepush = true)
