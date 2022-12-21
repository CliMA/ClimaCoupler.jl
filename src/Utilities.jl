#=
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
=#
module Utilities

using ClimaCore: ClimaCore, Fields, Spaces, Domains, Meshes, Topologies
using ClimaComms

export CoupledSimulation, float_type, heaviside, swap_space!, create_space

"""
Stores information needed to run a simulation with the coupler. 
"""
struct CoupledSimulation{FT, S, D, B, FV, P, E}
    tspan::S
    dates::D
    boundary_space::B
    fields::FV
    parsed_args::P
    conservation_checks::E
    t::FT
    Δt_cpl::FT
    surface_masks::NamedTuple
    model_sims::NamedTuple
    mode::NamedTuple
    monthly_3d_diags::NamedTuple
    monthly_2d_diags::NamedTuple
end

CoupledSimulation{FT}(args...) where {FT} = CoupledSimulation{FT, typeof.(args[1:6])...}(args...)
float_type(::CoupledSimulation{FT}) where {FT} = FT

"""
    heaviside(var)

Implements the heaviside step function multiplied by `var`, returning 0 for negative inputs
or the input value itself for non-negative inputs.

# Arguments
- `var`: [Integer or Float] value to apply heaviside to.
"""

heaviside(var) = var < 0 ? 0 : var

"""
    swap_space!(field::Fields.Field, new_space)

Remap the values of a field onto a new space.

# Arguments
- `field`: [Fields.Field] to be remapped to new space.
- `new_space`: [Spaces.Space] to remap `field` to.
"""
function swap_space!(field::Fields.Field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

end
