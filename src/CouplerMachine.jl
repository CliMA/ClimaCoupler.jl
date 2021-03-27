"""
    CouplerMachine

Primitive coupling module sufficient for initial atmos-ocean-land coupled simulation.
"""
module CouplerMachine

using ClimateMachine
using ClimateMachine.DGMethods

include("CplModel.jl")
include("CplState.jl")
include("CplSolver.jl")

end
