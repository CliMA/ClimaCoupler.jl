"""
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
"""
module Utilities

import ClimaComms
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
    field_out = Fields.Field(Fields.field_values(field_in), axes(field_out))
end

"""
    get_device(parsed_args) 

Returns the device on which the model is being run 

# Arguments 
- `parsed_args`: dictionary containing a "device" flag which decides which device to run on 
"""
function get_device(parsed_args)
    if parsed_args["device"] == "auto"
        return ClimaComms.device()
    elseif parsed_args["device"] == "CUDADevice"
        return ClimaComms.CUDADevice()
    elseif parsed_args["device"] == "CPUMultiThreaded" || Threads.nthreads() > 1
        return ClimaComms.CPUMultiThreaded()
    else
        return ClimaComms.CPUSingleThreaded()
    end
end


"""
    get_comms_context(parsed_args) 

Sets up the appropriate ClimaComms context for the device the model is to be run on 

# Arguments 
`parsed_args`: dictionary containing a "device" flag whcih decides which device context is needed
"""
function get_comms_context(parsed_args)
    device = get_device(parsed_args)
    comms_ctx = ClimaComms.context(device)
    ClimaComms.init(comms_ctx)

    @info "Running on $(nameof(typeof(device)))."

    if comms_ctx isa ClimaComms.SingletonCommsContext
        @info "Setting up single-process ClimaCoupler run"
    else
        @info "Setting up distributed ClimaCoupler run" nprocs = ClimaComms.nprocs(comms_ctx)
    end

    return comms_ctx
end

end # module
