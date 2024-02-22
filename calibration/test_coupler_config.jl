import CalibrateAtmos
using Test

include("coupler_interface.jl")

# Tests for ensuring CalibrateAtmos sets AtmosConfig correctly.

member_path = joinpath("test_output", "iteration_001", "member_001")
file_path = joinpath(member_path, "parameters.toml")
mkpath(dirname(file_path))
touch(file_path)

config_dict = Dict{Any, Any}(
    "dt_save_to_disk" => "100days",
    "output_dir" => "test_output",
)

coupler_config = get_coupler_sim(1, 1, "amip_coupled")
(; parsed_args) = coupler_config

@testset "Atmos Configuration" begin
    @test parsed_args["moist"] == "equil"
    @test parsed_args["toml"] == [file_path]
    @test parsed_args["output_dir"] == member_path
end

rm(file_path)
