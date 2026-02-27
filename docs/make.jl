using ClimaCoupler
using Documenter, Literate
# Import packages needed to load plotting extensions
import CairoMakie, ClimaCoreMakie, GeoMakie, Makie, Poppler_jll, Printf, Oceananigans
# Import packages for ClimaCouplerCMIPExt
import ClimaOcean, ClimaSeaIce, KernelAbstractions
# Import packages for ClimaCouplerClimaLandExt
import ClimaLand
# Import packages for ClimaCouplerClimaAtmosExt
import ClimaAtmos

const COUPLER_DIR = joinpath(@__DIR__, "..")
const EXPERIMENTS_DIR = joinpath(@__DIR__, "..", "experiments")
const OUTPUT_DIR = joinpath(@__DIR__, "src/generated")

# tutorials & experiments
# - generate tutorial files:

# sea breeze tutorial
TUTORIAL_DIR_SB = joinpath(EXPERIMENTS_DIR, "ClimaCore/sea_breeze/")
TUTORIAL_DIR_AMIP = joinpath(EXPERIMENTS_DIR, "ClimaEarth/")

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

# pages layout
example_pages = [
    "AMIP" => map(
        s -> "generated/amip/$(s)",
        readdir(joinpath(@__DIR__, "src/generated/amip")),
    ),
    "Sea Breeze" => map(
        s -> "generated/sea_breeze/$(s)",
        readdir(joinpath(@__DIR__, "src/generated/sea_breeze")),
    ),
]
interface_pages = [
    "interfacer.md",
    "input.md",
    "simcoordinator.md",
    "timemanager.md",
    "fieldexchanger.md",
    "fluxcalculator.md",
    "checkpointer.md",
    "conservation.md",
    "utilities.md",
]

output_pages = ["simoutput.md", "plotting.md", "leaderboard.md"]

developer_pages = ["contributing.md", "debugging.md"]

pages = Any[
    "Home" => "index.md",
    "Running a simulation" => "running.md",
    "The CoupledSimulation object" => "coupledsimulation.md",
    "Available component models" => "models.md",
    "Available simulation types" => "simtypes.md",
    "Coupler interface" => interface_pages,
    "Simulation output" => output_pages,
    "Examples" => example_pages,
    "Parameter calibration" => "calibrationtools.md",
    "Developer docs" => developer_pages,
]

makedocs(
    modules = [
        ClimaCoupler,
        Base.get_extension(ClimaCoupler, :ClimaCouplerMakieExt),
        Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPMakieExt),
        Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt),
        Base.get_extension(ClimaCoupler, :ClimaCouplerClimaLandExt),
        Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt),
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
