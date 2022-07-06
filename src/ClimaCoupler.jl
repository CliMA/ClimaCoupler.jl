"""
    ClimaCoupler

Coupling module sufficient for initial atmos-ocean-land coupled simulation.
"""
module ClimaCoupler

using ClimaCore
import ClimaCore: Fields, Operators
import Logging

# Disable errors from broadcasting with mismatched spaces
Fields.allow_mismatched_diagonalized_spaces() = true
Operators.allow_mismatched_fd_spaces() = true

function __init__()
    Logging.disable_logging(Logging.Warn) # disable warnings (for broadcasting with mismatched spaces)
end

include("CoupledSimulations/clock.jl")
include("CoupledSimulations/coupled_simulation.jl")
include("CouplerState/coupler_state.jl")

end
