#=
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
=#
module Utilities

using ClimaCore: Fields, Spaces

export CoupledSimulation, float_type_cs, swap_space!


"""
Stores information needed to run a simulation with the coupler. 
"""
struct CoupledSimulation{FT, X, S, D, B, FV, P, E}
    comms_ctx::X
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

CoupledSimulation{FT}(args...) where {FT} = CoupledSimulation{FT, typeof.(args[1:7])...}(args...)
float_type_cs(::CoupledSimulation{FT}) where {FT} = FT

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
