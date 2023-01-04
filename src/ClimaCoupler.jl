"""
    ClimaCoupler

Coupling module sufficient for initial atmos-ocean-land coupled simulation.
"""
module ClimaCoupler

include("CoupledSimulations/clock.jl")
include("CoupledSimulations/coupled_simulation.jl")
include("CouplerState/coupler_state.jl")
include("../test/TestHelper.jl")
include("Utilities.jl")
include("Regridder.jl")
include("ConservationChecker.jl")
include("CallbackManager.jl")
include("BCReader.jl")
include("Diagnostics.jl")


end
