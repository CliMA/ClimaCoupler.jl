import ArgParse
function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        # ClimaCoupler flags
        "--run_name"
        help = "Name of this run."
        arg_type = String
        "--dt_cpl"
        help = " Coupling time step in seconds"
        arg_type = Int
        default = 400
        "--anim"
        help = "Boolean flag indicating whether to make animations"
        arg_type = Bool
        default = false
        "--energy_check"
        help = "Boolean flag indicating whether to check energy conservation"
        arg_type = Bool
        default = false
        "--mode_name"
        help = "Mode of coupled simulation. [`amip`, `slabplanet`]"
        arg_type = String
        default = "amip"
        "--mono_surface"
        help = "Boolean flag indicating whether (1st order) monotone and conservative remapping is applied."
        arg_type = Bool
        default = false
        "--turb_flux_partition"
        help = "Method to partition turbulent fluxes. [`PartitionedStateFluxes`, `CombinedStateFluxes`]"
        arg_type = String
        default = "CombinedStateFluxes"
        "--monthly_checkpoint" # TODO generalize to any frequency
        help = "Boolean flag indicating whether to checkpoint monthly"
        arg_type = Bool
        default = false
        "--restart_dir"
        help = "Directory containing restart files"
        arg_type = String
        default = "unspecified"
        "--restart_t"
        help = "Restart time"
        arg_type = Int
        default = 0
        "--config_file"
        help = "A yaml file used to set the configuration of the coupled model"
        "--print_config_dict"
        help = "Boolean flag indicating whether to print the final configuration dictionary"
        arg_type = Bool
        default = true
        "--FLOAT_TYPE"
        help = "Floating point precision  [`Float64` (default), `Float32`]"
        arg_type = String
        default = "Float64"
        # ClimaAtmos specific
        "--surface_setup"
        help = "Triggers ClimaAtmos into the coupled mode [`PrescribedSurface` (default)]" # retained here for standalone Atmos benchmarks
        arg_type = String
        default = "PrescribedSurface"
        "--atmos_config_file"
        help = "A yaml file used to set the atmospheric model configuration. If nothing is specified, the default configuration is used."
        # ClimaLSM specific
        "--land_albedo_type"
        help = "Access land surface albedo information from data file. [`function`, `map_static`, `map_temporal`]"
        arg_type = String
        default = "map_static" # to be replaced by land config file, when available
        "--land_domain_type"
        help = "Type of land domain. [`sphere` (default), `single_column`]"
        arg_type = String
        default = "sphere"
    end
    return s
end

parse_commandline(s) = ArgParse.parse_args(ARGS, s)

function cli_defaults(s::ArgParse.ArgParseSettings)
    defaults = Dict()
    # TODO: Don't use ArgParse internals
    for arg in s.args_table.fields
        defaults[arg.dest_name] = arg.default
    end
    return defaults
end

"""
    job_id_from_parsed_args(
        s::ArgParseSettings,
        parsed_args = ArgParse.parse_args(ARGS, s)
    )
Returns a unique name (`String`) given
 - `s::ArgParse.ArgParseSettings` The arg parse settings
 - `parsed_args` The parse arguments
The `ArgParseSettings` are used for truncating
this string based on the default values.
"""
job_id_from_parsed_args(s, parsed_args = ArgParse.parse_args(ARGS, s)) =
    job_id_from_parsed_args(cli_defaults(s), parsed_args)

function job_id_from_parsed_args(defaults::Dict, parsed_args)
    _parsed_args = deepcopy(parsed_args)
    s = ""
    warn = false
    for k in keys(_parsed_args)
        # Skip defaults to alleviate verbose names
        !haskey(defaults, k) && continue
        defaults[k] == _parsed_args[k] && continue

        if _parsed_args[k] isa String
            # We don't need keys if the value is a string
            # (alleviate verbose names)
            s *= _parsed_args[k]
        elseif _parsed_args[k] isa Int
            s *= k * "_" * string(_parsed_args[k])
        elseif _parsed_args[k] isa AbstractFloat
            warn = true
        else
            s *= k * "_" * string(_parsed_args[k])
        end
        s *= "_"
    end
    s = replace(s, "/" => "_")
    s = strip(s, '_')
    warn && @warn "Truncated job ID:$s may not be unique due to use of Real"
    return s
end


"""
    print_repl_script(str::String)
Generate a block of code to run a particular
buildkite job given the `command:` string.
Example:
"""
function print_repl_script(str)
    ib = """"""
    ib *= """\n"""
    ib *= """using Revise; include("src/utils/cli_options.jl");\n"""
    ib *= """\n"""
    ib *= """parsed_args = parse_commandline(argparse_settings());\n"""
    parsed_args = parsed_args_from_command_line_flags(str)
    for (flag, val) in parsed_args
        if val isa AbstractString
            ib *= "parsed_args[\"$flag\"] = \"$val\";\n"
        else
            ib *= "parsed_args[\"$flag\"] = $val;\n"
        end
    end
    ib *= """\n"""
    ib *= """include("examples/hybrid/driver.jl")\n"""
    println(ib)
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

parsed_args_from_ARGS(ARGS, parsed_args = Dict()) = parsed_args_from_ARGS_string(strip(join(ARGS, " ")), parsed_args)

parsed_args_from_command_line_flags(str, parsed_args = Dict()) =
    parsed_args_from_ARGS_string(strip(last(split(str, ".jl"))), parsed_args)

function parsed_args_from_ARGS_string(str, parsed_args = Dict())
    str = replace(str, "    " => " ", "   " => " ", "  " => " ")
    parsed_args_list = split(str, " ")
    parsed_args_list == [""] && return parsed_args
    @assert iseven(length(parsed_args_list))
    parsed_arg_pairs = map(1:2:(length(parsed_args_list) - 1)) do i
        Pair(parsed_args_list[i], strip(parsed_args_list[i + 1], '\"'))
    end
    function parse_arg(val)
        for T in (Bool, Int, Float32, Float64)
            try
                return parse(T, val)
            catch
            end
        end
        return String(val) # string
    end
    for (flag, val) in parsed_arg_pairs
        parsed_args[replace(flag, "--" => "")] = parse_arg(val)
    end
    return parsed_args
end

"""
    parsed_args_per_job_id()
    parsed_args_per_job_id(buildkite_yaml)

A dict of `parsed_args` to run the ClimaAtmos driver
whose keys are the `job_id`s from buildkite yaml.

# Example

To run the `sphere_aquaplanet_rhoe_equilmoist_allsky`
buildkite job from the standard buildkite pipeline, use:
```
using Revise; include("examples/hybrid/cli_options.jl");
dict = parsed_args_per_job_id();
parsed_args = dict["sphere_aquaplanet_rhoe_equilmoist_allsky"];
include("examples/hybrid/driver.jl")
```
"""
function parsed_args_per_job_id(; trigger = "driver.jl")
    cc_dir = joinpath(@__DIR__, "..", "..", "..")
    buildkite_yaml = joinpath(cc_dir, ".buildkite", "pipeline.yml")
    parsed_args_per_job_id(buildkite_yaml; trigger)
end

function parsed_args_per_job_id(buildkite_yaml; trigger = "driver.jl")
    buildkite_commands = readlines(buildkite_yaml)
    filter!(x -> occursin(trigger, x), buildkite_commands)

    @assert length(buildkite_commands) > 0 # sanity check
    result = Dict()
    for bkcs in buildkite_commands
        default_parsed_args = parse_commandline(argparse_settings())
        job_id = first(split(last(split(bkcs, "--run_name ")), " "))
        job_id = strip(job_id, '\"')
        result[job_id] = parsed_args_from_command_line_flags(bkcs, default_parsed_args)
    end
    return result
end

function non_default_command_line_flags_parsed_args(parsed_args)
    default_parsed_args = parse_commandline(argparse_settings())
    s = ""
    for k in keys(parsed_args)
        default_parsed_args[k] == parsed_args[k] && continue
        s *= "--$k $(parsed_args[k]) "
    end
    return rstrip(s)
end
