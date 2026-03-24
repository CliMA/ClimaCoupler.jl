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
include("Models.jl")
include("CalibrationTools.jl")

# Import key functions from submodules to re-export at top level
import ..Interfacer: CoupledSimulation
import ..SimCoordinator: run!, step!, setup_and_run
import ..Plotting: postprocess

# Export all modules and key functions
export CalibrationTools,
    ConservationChecker,
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Input,
    Interfacer,
    Models,
    Plotting,
    SimCoordinator,
    SimOutput,
    TimeManager,
    Utilities,
    CoupledSimulation,
    run!,
    step!,
    setup_and_run,
    postprocess

end
