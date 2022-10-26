#=
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
=#
module Utilities

using ClimaCore: ClimaCore, Fields, Spaces, Domains, Meshes, Topologies
using ClimaComms

export CoupledSimulation, heaviside, swap_space!, create_space

"""
Stores information needed to run a simulation with the coupler. 
"""
struct CoupledSimulation{I, F, S, D, B, T, FV, P, E}
    Δt_cpl::I
    t::F
    tspan::S
    dates::D
    boundary_space::B
    FT::T
    surface_masks::NamedTuple
    fields::FV
    model_sims::NamedTuple
    mode::NamedTuple
    parsed_args::P
    monthly_3d_diags::NamedTuple
    monthly_2d_diags::NamedTuple
    conservation_checks::E
end

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
