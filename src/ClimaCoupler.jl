"""
    ClimaCoupler

Module for atmos-ocean-land coupled simulations.
"""
module ClimaCoupler

include("Interfacer.jl")
include("CoupledSimulations/clock.jl")
include("CoupledSimulations/coupled_simulation.jl")
include("CouplerState/coupler_state.jl")
include("../test/TestHelper.jl")
include("Utilities.jl")
include("TimeManager.jl")
include("Regridder.jl")
include("ConservationChecker.jl")
include("BCReader.jl")
include("FluxCalculator.jl")
include("FieldExchanger.jl")
include("Diagnostics.jl")
include("PostProcessor.jl")

end
