"""
    ClimaCoupler

Coupling module sufficient for initial atmos-ocean-land coupled simulation.
"""
module ClimaCoupler

using DocStringExtensions

@template DEFAULT =
    """
    (DEFAULT)
    $(TYPEDEF)
    $(DOCSTRING)
    """
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
    """

"test dummy struct"
mutable struct Bongo end

"test struct 2"
mutable struct Alpha{S} end



include("CoupledSimulations/clock.jl")
include("CoupledSimulations/coupled_simulation.jl")
include("CouplerState/coupler_state.jl")

end
