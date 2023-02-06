"""
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
"""
module Utilities

using ClimaCore: Fields, Spaces

export CoupledSimulation, float_type, swap_space!

"""
    AbstractSimulation

An abstract super-type representing a simulation.
"""
abstract type AbstractSimulation{FT} end

"""
    CoupledSimulation
Stores information needed to run a simulation with the coupler.
"""
struct CoupledSimulation{
    FT <: Real,
    X,
    D,
    B,
    FV,
    P,
    E,
    TS,
    TI <: Real,
    DTI <: Real,
    NTSM <: NamedTuple,
    NTMS <: NamedTuple,
    NTM <: NamedTuple,
}
    comms_ctx::X
    dates::D
    boundary_space::B
    fields::FV
    parsed_args::P
    conservation_checks::E
    tspan::TS
    t::TI
    Δt_cpl::DTI
    surface_masks::NTSM
    model_sims::NTMS
    mode::NTM
    diagnostics::Tuple
end

CoupledSimulation{FT}(args...) where {FT} = CoupledSimulation{FT, typeof.(args[1:12])...}(args...)

"""
    float_type(::CoupledSimulation)

Return the floating point type backing `T`: `T` can either be an object or a type.
"""
float_type(::CoupledSimulation{FT}) where {FT} = FT

"""
    swap_space!(field_out::Fields.Field, field_in::Fields.Field)

Remap the values of a field onto a new space.

# Arguments
- `field_in`: [Fields.Field] to be remapped to new space.
- `field_out`: [Fields.Field] to remap `field_in` to.
"""
function swap_space!(field_out, field_in::Fields.Field)
    parent(field_out) .= parent(field_in)
    return field_out
end

end # module
