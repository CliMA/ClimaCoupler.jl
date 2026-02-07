using ClimaCoupler
using Documenter, Literate
# Import packages needed to load plotting extensions
import CairoMakie, ClimaCoreMakie, GeoMakie, Makie, Poppler_jll, Printf, Oceananigans
# Import packages for ClimaCouplerCMIPExt
import ClimaOcean, ClimaSeaIce, KernelAbstractions
# Import packages for ClimaCouplerClimaLandExt
import ClimaLand, NCDatasets

const COUPLER_DIR = joinpath(@__DIR__, "..")
const EXPERIMENTS_DIR = joinpath(@__DIR__, "..", "experiments")
const OUTPUT_DIR = joinpath(@__DIR__, "src/generated")

# tutorials & experiments
# - generate tutorial files:

# sea breeze tutorial
TUTORIAL_DIR_SB = joinpath(EXPERIMENTS_DIR, "ClimaCore/sea_breeze/")
TUTORIAL_DIR_AMIP = joinpath(EXPERIMENTS_DIR, "ClimaEarth/")

# Pkg.activate(TUTORIAL_DIR)
# Pkg.instantiate()
# include(joinpath(TUTORIAL_DIR, "run.jl"))
# Literate.markdown(joinpath(TUTORIAL_DIR, tutorial_name), OUTPUT_DIR; execute = true, documenter = false)

# execute Literate on all julia files
tutorial_files_sb = filter(x -> last(x, 3) == ".jl", readdir(TUTORIAL_DIR_SB))
tutorial_files_amip = filter(x -> last(x, 11) == "run_amip.jl", readdir(TUTORIAL_DIR_AMIP))

# Literate generates markdown files and stores them in docs/src/generated/sea_breeze
map(
    x -> Literate.markdown(
        joinpath(TUTORIAL_DIR_SB, x),
        joinpath(OUTPUT_DIR, "sea_breeze");
        execute = false,
        documenter = false,
    ),
    tutorial_files_sb,
)

map(
    x -> Literate.markdown(
        joinpath(TUTORIAL_DIR_AMIP, x),
        joinpath(OUTPUT_DIR, "amip");
        execute = false,
        documenter = false,
    ),
    tutorial_files_amip,
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
experiment_pages = [
    "Sea Breeze" => map(
        s -> "generated/sea_breeze/$(s)",
        readdir(joinpath(@__DIR__, "src/generated/sea_breeze")),
    ),
    "AMIP" => map(
        s -> "generated/amip/$(s)",
        readdir(joinpath(@__DIR__, "src/generated/amip")),
    ),
]
interface_pages = [
    "input.md",
    "checkpointer.md",
    "conservation.md",
    "fieldexchanger.md",
    "fluxcalculator.md",
    "interfacer.md",
    "models.md",
    "simcoordinator.md",
    "timemanager.md",
    "utilities.md",
    "simoutput.md",
    "plotting.md",
    "calibrationtools.md",
]
performance_pages = ["performance.md"]

output_pages = ["diagnostics.md", "leaderboard.md"]

pages = Any[
    "Home" => "index.md",
    "Examples" => experiment_pages,
    "Coupler Interface" => interface_pages,
    "Performance" => performance_pages,
    "Model Output" => output_pages,
    "Contributing" => "contributing.md",
]


makedocs(
    modules = [
        ClimaCoupler,
        Base.get_extension(ClimaCoupler, :ClimaCouplerMakieExt),
        Base.get_extension(ClimaCoupler, :ClimaCouplerOceananigansMakieExt),
        Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt),
        Base.get_extension(ClimaCoupler, :ClimaCouplerClimaLandExt),
    ],
    authors = "Climate Modelling Alliance",
    sitename = "ClimaCoupler.jl",
    format = Documenter.HTML(),
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/ClimaCoupler.jl.git",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
    versions = ["stable" => "v^", "v#.#.#", "dev" => "main"],
)
