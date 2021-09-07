"""
    CouplerMachine

Coupling module sufficient for initial atmos-ocean-land coupled simulation.
"""
module CouplerMachine

include("CoupledSimulations/coupled_simulation.jl")
include("CouplerState/coupler_state.jl")

end
