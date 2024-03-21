import DelimitedFiles as DLM
using Statistics
import ClimaCoupler
import ArgParse

function argparse_settings()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! s begin
        "--cpu_run_name"
        help = "The name of the CPU run we want to compare. User must specify CPU and/or GPU run name."
        arg_type = String
        default = nothing
        "--gpu_run_name"
        help = "The name of the GPU run we want to compare."
        arg_type = String
        default = nothing
        "--mode_name"
        help = "The mode of the simulations being compared (`slabplanet` or `AMIP`)."
        arg_type = String
        default = nothing
        "--coupler_output_dir"
        help = "Directory to save output files."
        arg_type = String
        default = "experiments/AMIP/output"
    end
    return s
end

# Read in CPU and GPU run name info from command line
parsed_args = parse_commandline(argparse_settings())
cpu_run_name = parsed_args["cpu_run_name"]
gpu_run_name = parsed_args["gpu_run_name"]
if isnothing(cpu_run_name) && isnothing(gpu_run_name)
    error("Must pass CPU and/or GPU run name to compare them.")
elseif isnothing(gpu_run_name)
    gpu_run_name = "gpu_" * cpu_run_name
elseif isnothing(cpu_run_name)
    cpu_run_name = gpu_run_name[5:end]
end

# Read in mode name from command line (or retrieve from run name)
mode_name = parsed_args["mode_name"]
if isnothing(mode_name)
    mode_name =
        occursin("amip", cpu_run_name) ? "amip" :
        (occursin("slabplanet", cpu_run_name) ? "slabplanet" : error("Please provide a valid `mode_name`."))
end

# Construct CPU and GPU artifacts directories
output_dir = parsed_args["coupler_output_dir"]
cpu_artifacts_dir = joinpath(output_dir, joinpath(mode_name, cpu_run_name)) * "_artifacts"
gpu_artifacts_dir = joinpath(output_dir, joinpath(mode_name, gpu_run_name)) * "_artifacts"


# Read in and compare atmos state variables
cpu_atmos_state_center = DLM.readdlm(joinpath(cpu_artifacts_dir, "atmos_state_center_tend_cpu.txt"), ',')
gpu_atmos_state_center = DLM.readdlm(joinpath(gpu_artifacts_dir, "atmos_state_center_tend_gpu.txt"), ',')

@show abs(maximum(cpu_atmos_state_center .- gpu_atmos_state_center))
@show abs(median(cpu_atmos_state_center .- gpu_atmos_state_center))
@show abs(mean(cpu_atmos_state_center .- gpu_atmos_state_center))
@assert isapprox(cpu_atmos_state_center, gpu_atmos_state_center)

cpu_atmos_state_face = DLM.readdlm(joinpath(cpu_artifacts_dir, "atmos_state_face_tend_cpu.txt"), ',')
gpu_atmos_state_face = DLM.readdlm(joinpath(gpu_artifacts_dir, "atmos_state_face_tend_gpu.txt"), ',')

@show abs(maximum(cpu_atmos_state_face .- gpu_atmos_state_face))
@show abs(median(cpu_atmos_state_face .- gpu_atmos_state_face))
@show abs(mean(cpu_atmos_state_face .- gpu_atmos_state_face))
@assert isapprox(cpu_atmos_state_face, gpu_atmos_state_face)


# Read in and compare land state variables
cpu_land_state_3d = DLM.readdlm(joinpath(cpu_artifacts_dir, "land_state_3d_tend_cpu.txt"), ',')
gpu_land_state_3d = DLM.readdlm(joinpath(gpu_artifacts_dir, "land_state_3d_tend_gpu.txt"), ',')

@show abs(maximum(cpu_land_state_3d .- gpu_land_state_3d))
@show abs(median(cpu_land_state_3d .- gpu_land_state_3d))
@show abs(mean(cpu_land_state_3d .- gpu_land_state_3d))
@assert isapprox(cpu_land_state_3d, gpu_land_state_3d)

cpu_land_state_2d = DLM.readdlm(joinpath(cpu_artifacts_dir, "land_state_2d_tend_cpu.txt"), ',')
gpu_land_state_2d = DLM.readdlm(joinpath(gpu_artifacts_dir, "land_state_2d_tend_gpu.txt"), ',')

@show abs(maximum(cpu_land_state_2d .- gpu_land_state_2d))
@show abs(median(cpu_land_state_2d .- gpu_land_state_2d))
@show abs(mean(cpu_land_state_2d .- gpu_land_state_2d))
@assert isapprox(cpu_land_state_2d, gpu_land_state_2d)


# Read in ocean state variables (if not AMIP)
if !(mode_name == "amip")
    cpu_land_state = DLM.readdlm(joinpath(cpu_artifacts_dir, "land_state_tend_cpu.txt"), ',')
    gpu_land_state = DLM.readdlm(joinpath(gpu_artifacts_dir, "land_state_tend_gpu.txt"), ',')

    @show abs(maximum(cpu_ocean_state .- gpu_ocean_state))
    @show abs(median(cpu_ocean_state .- gpu_ocean_state))
    @show abs(mean(cpu_ocean_state .- gpu_ocean_state))
    @assert isapprox(cpu_ocean_state, gpu_ocean_state)
end
