"""
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
"""
module Utilities

import Artifacts
import ClimaComms
import ClimaCore as CC
import Logging

export swap_space!

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

    # Ensure that logging only happens on the root process when using multiple processes
    ClimaComms.iamroot(comms_ctx) || Logging.disable_logging(Logging.AboveMaxLevel)

    if comms_ctx isa ClimaComms.SingletonCommsContext
        @info "Setting up single-process ClimaCoupler run on device: $(nameof(typeof(device)))."
    else
        @info "Setting up distributed ClimaCoupler run on " nprocs = ClimaComms.nprocs(comms_ctx) device = "$(nameof(typeof(device)))"
    end

    return comms_ctx
end

"""
    show_memory_usage()

Display and return the maximum resident set size (RSS) memory footprint on the
CPU of this process since it began.
"""
function show_memory_usage()
    cpu_max_rss_GB = ""
    cpu_max_rss_GB = string(round(Sys.maxrss() / 1e9, digits = 3)) * " GiB"
    @info cpu_max_rss_GB
    return cpu_max_rss_GB
end

end # module
