# helpers with boiler plate code for IO operations, useful for all ClimaEarth drivers.
import ClimaUtilities.OutputPathGenerator: generate_output_path


"""
    setup_output_dirs(; output_dir = nothing, artifacts_dir = nothing, comms_ctx)

Create output directories for the experiment. If `comms_ctx` is provided, only the root process will create the directories.
By default, the regrid directory is created as a temporary directory inside the output directory,
and the artifacts directory is created inside the output directory with the suffix `_artifacts`.

`ClimaUtilities.OutputPathGenerator` is used so that simulations can be re-run and re-started.
The output path looks like:
```
output_dir/
├── output_0000/
│   └── ... component model outputs in their folders ...
├── output_0001/
│   └── ... component model outputs in their folders ...
├── output_0002/
│   └── ... component model outputs in their folders ...
└── output_active -> output_0002/
```

# Arguments
- `output_dir::String`: The directory where the output files will be stored. Default is the current directory.
- `regrid_dir::String`: The directory where the regridded files will be stored. Default is `output_dir/regrid_tmp/`.
- `artifacts_dir::String`: The directory where the artifacts will be stored. Default is `output_dir_artifacts/`.
- `comms_ctx::Union{Nothing, ClimaComms.AbstractCommsContext}`: The communicator context. If provided, only the root process will create the directories.

# Returns
- A tuple with the paths to the output, regrid, and artifacts directories.
"""
function setup_output_dirs(; output_dir = pwd(), artifacts_dir = joinpath(output_dir, "artifacts"), comms_ctx)
    output_dir = generate_output_path(output_dir)

    @info(output_dir)
    regrid_dir = nothing
    if ClimaComms.iamroot(comms_ctx)
        mkpath(artifacts_dir)
        regrid_dir = mktempdir(output_dir, prefix = "regrid_tmp_")
    end
    regrid_dir = ClimaComms.bcast(comms_ctx, regrid_dir)

    return (; output = output_dir, artifacts = artifacts_dir, regrid = regrid_dir)
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

"""
    get_period(t_start, t_end)

Return a period that depends on the duration of the simulation. Used for diagnostics.
"""
function get_period(t_start, t_end)
    sim_duration = t_end - t_start
    secs_per_day = 86400
    if sim_duration >= 90 * secs_per_day
        # if duration >= 90 days, take monthly means
        period = "1months"
        calendar_dt = Dates.Month(1)
    elseif sim_duration >= 30 * secs_per_day
        # if duration >= 30 days, take means over 10 days
        period = "10days"
        calendar_dt = Dates.Day(10)
    elseif sim_duration >= secs_per_day
        # if duration >= 1 day, take daily means
        period = "1days"
        calendar_dt = Dates.Day(1)
    else
        # if duration < 1 day, take hourly means
        period = "1hours"
        calendar_dt = Dates.Hour(1)
    end
    return (period, calendar_dt)
end
