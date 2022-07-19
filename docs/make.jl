using ClimaCoupler
using Documenter, Literate, Pkg

const COUPLER_DIR = joinpath(@__DIR__, "..")
const EXPERIMENTS_DIR = joinpath(@__DIR__, "..", "experiments")
const OUTPUT_DIR = joinpath(@__DIR__, "src/generated")

# tutorials & experiments
# - generate tutorial files:

# sea breeze tutorial
TUTORIAL_DIR = joinpath(EXPERIMENTS_DIR, "ClimaCore/sea_breeze/")
# Pkg.activate(TUTORIAL_DIR)
# Pkg.instantiate()
# include(joinpath(TUTORIAL_DIR, "run.jl"))
# Literate.markdown(joinpath(TUTORIAL_DIR, tutorial_name), OUTPUT_DIR; execute = true, documenter = false)
# execute Literate on all julia files
tutorial_files = filter(x -> last(x, 3) == ".jl", readdir(TUTORIAL_DIR))
# Literate generates markdown files and stores them in docs/src/generated/sea_breeze
map(
    x -> Literate.markdown(
        joinpath(TUTORIAL_DIR, x),
        joinpath(OUTPUT_DIR, "sea_breeze");
        execute = false,
        documenter = false,
    ),
    tutorial_files,
)

# - move tutorial files to docs/src
# IMAGE_DIR = joinpath(TUTORIAL_DIR, "images/")
# files = readdir(IMAGE_DIR)
# png_files = filter(endswith(".png"), files)
# for file in png_files
#     mkpath(joinpath(OUTPUT_DIR, "images/"))
#     cp(joinpath(IMAGE_DIR, file), joinpath(OUTPUT_DIR, "images/", file), force = true)
# end

# pages layout
experiment_pages =
    ["Sea Breeze" => map(s -> "generated/sea_breeze/$(s)", readdir(joinpath(@__DIR__, "src/generated/sea_breeze")))]
interface_pages = ["couplerstate.md", "timestepping.md"]

pages = Any["Home" => "index.md", "Examples" => experiment_pages, "Coupler Interface" => interface_pages]


makedocs(sitename = "ClimaCoupler.jl", format = Documenter.HTML(), modules = [ClimaCoupler], pages = pages)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# deploydocs(repo = "<github.com/CliMA/ClimaCoupler.jl.git>", push_preview = true, devbranch = "main", forcepush = true)
