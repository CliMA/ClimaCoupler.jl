"""
    CouplerMachine

Coupling module sufficient for initial atmos-ocean-land coupled simulation.
"""
module CouplerMachine

include("CplSimulations/CoupledSimulation.jl")
include("CplSimulations/AOLSimulation.jl")
include("CplState/CouplerState.jl")

end
