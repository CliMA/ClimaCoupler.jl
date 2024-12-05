"""
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
"""
module Utilities

import ClimaComms
import ClimaCore as CC
import Logging
import ClimaUtilities.OutputPathGenerator: generate_output_path

export binary_mask, swap_space!, get_device, get_comms_context, show_memory_usage, setup_output_dirs, time_to_seconds

"""
    binary_mask(var, threshold)

Converts a number `var` to 1, if `var` is greater or equal than a given `threshold` value,
or 0 otherwise, keeping the same type.

# Arguments
- `var`: [FT] value to be converted.
- `threshold`: [FT] cutoff value for conversions.
"""
binary_mask(var, threshold) = var >= threshold ? one(var) : zero(var)

"""
    binary_mask(var)

Converts a number `var` to 1, if `var` is greater or equal than `eps(FT)`,
or 0 otherwise, keeping the same type.

# Arguments
- `var`: [FT] value to be converted.
"""
binary_mask(var) = binary_mask(var, eps(eltype(var)))

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
    get_device(config_dict)

Returns the device on which the model is being run

# Arguments
- `config_dict`: dictionary containing a "device" flag which decides which device to run on
"""
function get_device(config_dict)
    if config_dict["device"] == "auto"
        return ClimaComms.device()
    elseif config_dict["device"] == "CUDADevice"
        return ClimaComms.CUDADevice()
    elseif config_dict["device"] == "CPUMultiThreaded" || Threads.nthreads() > 1
        return ClimaComms.CPUMultiThreaded()
    else
        return ClimaComms.CPUSingleThreaded()
    end
end


"""
    get_comms_context(config_dict)

Sets up the appropriate ClimaComms context for the device the model is to be run on,
choosing from the following options:
    - CPU single threaded
    - CPU with MPI
    - GPU

If no device is passed to `ClimaComms.context()` then `ClimaComms` automatically
selects the device from which this code is called.

# Arguments
`config_dict`: dictionary containing a "device" flag whcih decides which device context is needed
"""
function get_comms_context(config_dict)
    device = get_device(config_dict)
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

"""
    setup_output_dirs(; output_dir = nothing, artifacts_dir = nothing, comms_ctx)

Create output directories for the experiment. If `comms_ctx` is provided, only the root process will create the directories.
By default, the regrid directory is created as a temporary directory inside the output directory,
and the artifacts directory is created inside the output directory with the name `artifacts/`.

`ClimaUtilities.OutputPathGenerator` is used so that simulations can be re-run and re-started.
The output path looks like:
```
coupler_output_dir_amip/
├── checkpoints
│       └── checkpoints for the various models
├── artifacts
│       └── plots produced by the postporcessing step
├── output_0000/
│   ├── atmos/
│   │   └── output of the atmos model
│   └── ocean/
│       └── output of the ocean model
├── output_0001/
│   └── ... component model outputs in their folders ...
├── output_0002/
│   └── ... component model outputs in their folders ...
└── output_active -> output_0002/
```

# Arguments
- `output_dir::String`: The directory where the output files will be stored. Default is the current directory.
- `regrid_dir::String`: The directory where the regridded files will be stored. Default is `output_dir/regrid_tmp/`.
- `checkpoint_dir::String`: The directory where the checkpoint files will be stored. Default is `output_dir/checkpoints/`.
- `artifacts_dir::String`: The directory where the artifacts will be stored. Default is `output_dir/artifacts/`.
- `comms_ctx::Union{Nothing, ClimaComms.AbstractCommsContext}`: The communicator context. If provided, only the root process will create the directories.

# Returns
- A tuple with the paths to the output, regrid, and artifacts directories.
"""
function setup_output_dirs(;
    output_dir = pwd(),
    artifacts_dir = joinpath(output_dir, "artifacts"),
    checkpoints_dir = joinpath(output_dir, "checkpoints"),
    comms_ctx,
)
    output_dir = generate_output_path(output_dir, context = comms_ctx)
    regrid_dir = nothing
    if ClimaComms.iamroot(comms_ctx)
        mkpath(artifacts_dir)
        mkpath(checkpoints_dir)
        regrid_dir = mktempdir(output_dir, prefix = "regrid_tmp_")
    end
    regrid_dir = ClimaComms.bcast(comms_ctx, regrid_dir)

    return (; output = output_dir, artifacts = artifacts_dir, regrid = regrid_dir, checkpoints = checkpoints_dir)
end

"""
    time_to_seconds(s::String)

Convert a string to seconds. The string should be in the format `numberunit`, where `unit` is one of `secs`, `mins`, `hours`, or `days`.

# Arguments
- `s::String`: The string to convert to seconds.

# Returns
- The number of seconds represented by the string.
"""
function time_to_seconds(s::String)
    factor = Dict("secs" => 1, "mins" => 60, "hours" => 60 * 60, "days" => 60 * 60 * 24)
    s == "Inf" && return Inf
    if count(occursin.(keys(factor), Ref(s))) != 1
        error("Bad format for flag $s. Examples: [`10secs`, `20mins`, `30hours`, `40days`]")
    end
    for match in keys(factor)
        occursin(match, s) || continue
        return parse(Float64, first(split(s, match))) * factor[match]
    end
    error("Uncaught case in computing time from given string.")
end


end # module
