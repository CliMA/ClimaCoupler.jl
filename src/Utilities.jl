"""
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
"""
module Utilities

using ClimaCore: Fields, Spaces

export swap_space!

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
