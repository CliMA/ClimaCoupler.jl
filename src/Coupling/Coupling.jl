"""
    Coupling

Primitive coupling module sufficient for initial atmos-ocean-land coupled simulation.
"""
module Coupling

export CplModel, CplState, put!, get, register_cpl_field!

using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.ODESolvers
using ClimateMachine.ODESolvers: AbstractODESolver
using ClimateMachine.Mesh.Grids: VerticalDirection

import ClimateMachine.Ocean.Domains: DiscontinuousSpectralElementGrid

include("CplModel.jl")
include("CplState.jl")

end
