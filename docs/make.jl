using ClimaCoupler
using Documenter, Literate
# Import packages needed to load plotting extensions
import CairoMakie, ClimaCoreMakie, GeoMakie, Makie, Poppler_jll, Printf, Oceananigans
# Import packages for ClimaCouplerCMIPExt
import ClimaOcean, ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
# Import packages for ClimaCouplerClimaLandExt
import ClimaLand
# Import packages for ClimaCouplerClimaAtmosExt
import ClimaAtmos

const COUPLER_DIR = joinpath(@__DIR__, "..")
const EXPERIMENTS_DIR = joinpath(@__DIR__, "..", "experiments")
const OUTPUT_DIR = joinpath(@__DIR__, "src/generated")

# tutorials & experiments
# - generate tutorial files:
TUTORIAL_DIR_AMIP = joinpath(EXPERIMENTS_DIR, "AMIP")

# execute Literate on all julia files
tutorial_files_amip = filter(x -> x == "run_simulation.jl", readdir(TUTORIAL_DIR_AMIP))

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

pages = Any[
    "Home" => "index.md",
    "Running a simulation" => "running.md",
    "The CoupledSimulation object" => "coupledsimulation.md",
    "Available component models" => "models.md",
    "Available simulation types" => "simtypes.md",
    "Coupler interface" => interface_pages,
    "Simulation output" => output_pages,
    "Examples" => example_pages,
    "Performance tips" => "performance.md",
    "Debugging tips" => "debugging.md",
    "Parameter calibration" => "calibrationtools.md",
    "Contributing to ClimaCoupler" => "contributing.md",
]


extensions = [
    Base.get_extension(ClimaCoupler, :ClimaCouplerMakieExt),
    Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPMakieExt),
    Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt),
    Base.get_extension(ClimaCoupler, :ClimaCouplerClimaLandExt),
    Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt),
]
@assert all(!isnothing, extensions) "Some extensions failed to load: $extensions"

makedocs(
    modules = [ClimaCoupler, extensions...],
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
