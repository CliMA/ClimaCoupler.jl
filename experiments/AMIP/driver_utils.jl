# useful functions to abstract the driver code

function setup_output_dirs(; output_dir = nothing, regrid_dir = nothing, artifacts_dir = nothing, comms_ctx = nothing)
    if output_dir === nothing
        output_dir = "."
    end
    if regrid_dir === nothing
        regrid_dir = joinpath(output_dir, "regrid_tmp/")
    end
    if artifacts_dir === nothing
        artifacts_dir = output_dir * "_artifacts"
    end

    if !isnothing(comms_ctx) && ClimaComms.iamroot(comms_ctx)
        @info(output_dir)
        mkpath(output_dir)
        mkpath(regrid_dir)
        mkpath(artifacts_dir)
    end

    !isnothing(comms_ctx) ? ClimaComms.barrier(comms_ctx) : nothing

    return (; output = output_dir, artifacts = artifacts_dir, regrid = regrid_dir)

end

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