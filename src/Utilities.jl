#=
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
=#
module Utilities

using ClimaCore: ClimaCore, Fields, Spaces

export CouplerSimulation, heaviside, swap_space!

"""
Stores information needed to run a simulation with the coupler. 
"""
struct CouplerSimulation{I, F, S, D, B, T, P}
    Δt_cpl::I
    t::F
    tspan::S
    dates::D
    boundary_space::B
    FT::T
    surface_masks::NamedTuple
    fields::NamedTuple
    model_sims::NamedTuple
    mode::NamedTuple
    parsed_args::P
    monthly_3d_diags::NamedTuple
    monthly_2d_diags::NamedTuple
end

"""
    heaviside(var)

Implements the heaviside step function, returning 0 for negative inputs
and the input value itself for non-negative inputs.

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
