"""
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
"""
module Utilities

using ClimaCore: Fields, Spaces
using Dates

export CoupledSimulation, float_type, swap_space!, current_date

"""
    swap_space!(space_out::CC.Spaces.AbstractSpace, field_in::CC.Fields.Field)

Remap the values of a field onto a new space.

# Arguments
- `space_out`: [CC.Spaces.AbstractSpace] The axes of the space we want to remap onto
- `field_in`: [CC.Fields.Field] to be remapped to new space.
"""
function swap_space!(space_out::CC.Spaces.AbstractSpace, field_in::CC.Fields.Field)
    field_out = CC.Fields.Field(CC.Fields.field_values(field_in), space_out)
    return field_out
end

"""
    current_date(cs::CoupledSimulation, t::Int)

Return the model date at the current timestep.

# Arguments
- `cs`: [CoupledSimulation] containing info about the simulation
- `t`: [Real] number of seconds since simulation began
"""
current_date(cs::CoupledSimulation, t::Real) = cs.dates.date0[1] + Dates.Second(t)

end # module
