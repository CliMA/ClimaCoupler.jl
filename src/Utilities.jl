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

export swap_space!, get_device, get_comms_context, show_memory_usage, setup_output_dirs, time_to_seconds

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
`config_dict`: dictionary containing a "device" flag which decides which device context is needed
"""
function get_comms_context(config_dict)
    device = get_device(config_dict)
    comms_ctx = ClimaComms.context(device)
    ClimaComms.init(comms_ctx)

    if pkgversion(ClimaComms) < v"0.6.6"
        # For older versions of ClimaComms, we have to manually ensure that logging only
        # happens on the root process when using multiple processes
        ClimaComms.iamroot(comms_ctx) || Logging.disable_logging(Logging.AboveMaxLevel)
    end

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
    @info "Memory in use: $(cpu_max_rss_GB)"
    return cpu_max_rss_GB
end

"""
    setup_output_dirs(output_dir = pwd(),
        artifacts_dir = joinpath(output_dir, "artifacts"),
        checkpoints_dir = joinpath(output_dir, "checkpoints"),
        regrid_dir = nothing,
        comms_ctx,
    )

Create output directories for the experiment. If `comms_ctx` is provided,
only the root process will create the directories.
By default, the artifacts and checkpoints directories are created inside the output
directory with the names `artifacts/` and `checkpoints/`.
The regrid directory is by default created as a temporary directory inside the output
directory and is automatically deleted when the process exits.


`ClimaUtilities.OutputPathGenerator` is used so that simulations can be re-run and re-started.
The output path looks like:
```
coupler_output_dir_amip/
├── checkpoints
│       └── checkpoints for the various models
├── artifacts
│       └── plots produced by the postporcessing step
├── regrid_tmp_<random_tempdir>/
│       └── temporary files used for regridding
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
- `regrid_dir::String`: The directory where the regridded files will be stored. Default is `output_dir/regrid_tmp_<random_tempdir>/`.
- `checkpoint_dir::String`: The directory where the checkpoint files will be stored. Default is `output_dir/checkpoints/`.
- `artifacts_dir::String`: The directory where the artifacts will be stored. Default is `output_dir/artifacts/`.
- `comms_ctx::Union{Nothing, ClimaComms.AbstractCommsContext}`: The communicator context. If provided, only the root process will create the directories.

# Returns
- A tuple with the paths to the output, artifacts, regrid, and checkpoints directories.
"""
function setup_output_dirs(;
    output_dir = pwd(),
    artifacts_dir = joinpath(output_dir, "artifacts"),
    checkpoints_dir = joinpath(output_dir, "checkpoints"),
    regrid_dir = nothing,
    comms_ctx,
)
    output_dir = generate_output_path(output_dir, context = comms_ctx)
    if ClimaComms.iamroot(comms_ctx)
        mkpath(artifacts_dir)
        mkpath(checkpoints_dir)
        # If no regrid_dir is provided, create a temporary directory
        regrid_dir = isnothing(regrid_dir) ? mktempdir(output_dir, prefix = "regrid_tmp_") : mkpath(regrid_dir)
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
