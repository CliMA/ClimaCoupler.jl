# helpers with boiler plate code for IO operations, useful for all ClimaEarth drivers.

"""
    setup_output_dirs(; output_dir = nothing, regrid_dir = nothing, artifacts_dir = nothing, comms_ctx)

Create output directories for the experiment. If `comms_ctx` is provided, only the root process will create the directories.

# Arguments
- `output_dir::String`: The directory where the output files will be stored. Default is the current directory.
- `regrid_dir::String`: The directory where the regridded files will be stored. Default is `output_dir/regrid_tmp/`.
- `artifacts_dir::String`: The directory where the artifacts will be stored. Default is `output_dir_artifacts/`.
- `comms_ctx::Union{Nothing, ClimaComms.AbstractCommsContext}`: The communicator context. If provided, only the root process will create the directories.

# Returns
- A tuple with the paths to the output, regrid, and artifacts directories.
"""
function setup_output_dirs(; output_dir = nothing, regrid_dir = nothing, artifacts_dir = nothing, comms_ctx)
    if output_dir === nothing
        output_dir = "."
    end
    if regrid_dir === nothing
        regrid_dir = joinpath(output_dir, "regrid_tmp/")
    end
    if artifacts_dir === nothing
        artifacts_dir = output_dir * "_artifacts"
    end

    @info(output_dir)
    if ClimaComms.iamroot(comms_ctx)
        mkpath(output_dir)
        mkpath(regrid_dir)
        mkpath(artifacts_dir)
    end

    ClimaComms.barrier(comms_ctx)

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
