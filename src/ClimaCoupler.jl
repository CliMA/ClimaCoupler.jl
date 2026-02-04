"""
    ClimaCoupler

Module for atmos-ocean-land coupled simulations.
"""
module ClimaCoupler

include("Interfacer.jl")
include("Utilities.jl")
include("TimeManager.jl")
include("ConservationChecker.jl")
include("FluxCalculator.jl")
include("FieldExchanger.jl")
include("Checkpointer.jl")
include("Input.jl")
include("SimOutput/SimOutput.jl")
include("Plotting.jl")
include("SimCoordinator.jl")

# Import run! and step! from SimCoordinator to re-export at top level
import ..SimCoordinator: run!, step!
# Import postprocess from Plotting to re-export at top level
import ..Plotting: postprocess

# Export all modules and key functions
export ConservationChecker,
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Input,
    Interfacer,
    Plotting,
    SimCoordinator,
    SimOutput,
    TimeManager,
    Utilities,
    run!,
    step!,
    postprocess

end
