"""
    ClimaCoupler

Module for atmos-ocean-land-ice coupled simulations.
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
include("Postprocessor/Postprocessor.jl")
include("SimCoordinator.jl")

# Re-export commonly used functions for convenience
export run!, step!  # from SimCoordinator
export postprocess  # from Postprocessor

end
