"""
    ClimaCoupler

Coupling module sufficient for initial atmos-ocean-land coupled simulation.
"""
module ClimaCoupler

using DocStringExtensions

@template FUNCTIONS =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    $(METHODLIST)
    """
@template METHODS =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    """
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    Fields:
    $(FIELDS)
    Contructors:
    $(METHODLIST)
    """

include("CoupledSimulations/clock.jl")
include("CoupledSimulations/coupled_simulation.jl")
include("CouplerState/coupler_state.jl")

end
