#=
    Unit tests for ClimaCoupler Utilities module
=#
import Test: @testset, @test, @test_throws
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCoupler: Utilities
import ClimaCore as CC
import Dates

# Initialize MPI context, in case
ClimaComms.init(ClimaComms.context())

@testset "parse_date / format_start_date" begin
    @test Utilities.parse_date("20000506") == Dates.DateTime(2000, 5, 6)
    @test Utilities.parse_date("20000506-0000") == Dates.DateTime(2000, 5, 6, 0, 0)
    @test Utilities.parse_date("20191231-1200") == Dates.DateTime(2019, 12, 31, 12, 0)
    @test Utilities.parse_date(Dates.DateTime(2019, 12, 31, 12)) ==
          Dates.DateTime(2019, 12, 31, 12)
    @test_throws ErrorException Utilities.parse_date("20000506-00000")
    @test_throws ErrorException Utilities.parse_date("")

    @test Utilities.format_start_date(Dates.DateTime(2019, 12, 31)) == "20191231"
    @test Utilities.format_start_date(Dates.DateTime(2019, 12, 31, 12)) == "20191231-1200"
end

for FT in (Float32, Float64)
    @testset "test comms_ctx" begin
        parsed_args = Dict("device" => "auto")

        this_device = ClimaComms.device()
        this_ctx = ClimaComms.context(this_device)

        @test typeof(Utilities.get_comms_context(parsed_args)) == typeof(this_ctx)

        if typeof(this_device) == typeof(ClimaComms.CUDADevice())
            parsed_args["device"] = "CUDADevice"
        elseif typeof(this_device) == typeof(ClimaComms.CPUMultiThreaded())
            parsed_args["device"] = "CPUMultiThreaded"
        else
            parsed_args["device"] = "CPUSingleThreaded"
        end

        ClimaComms.init(this_ctx)
        @test typeof(Utilities.get_comms_context(parsed_args)) == typeof(this_ctx)

        # Additional test calls for code coverage since Github Actions only exercises the SingleThreaded calculates
        # More meaningful testing performed on buildkite
        # Cannot test CUDADevice on CPU
        # Test other devices:
        parsed_args["device"] = "CPUMultiThreaded"
        @test typeof(Utilities.get_comms_context(parsed_args)) ==
              typeof(ClimaComms.context(ClimaComms.CPUMultiThreaded()))

        parsed_args["device"] = ""
        @test typeof(Utilities.get_comms_context(parsed_args)) ==
              typeof(ClimaComms.context(ClimaComms.CPUSingleThreaded()))
    end

    @testset "integral" begin
        space2d = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6), # in meters
            n_quad_points = 4,
            h_elem = 4,
        )
        ones2d = ones(space2d)

        space3d = CC.CommonSpaces.ExtrudedCubedSphereSpace(
            FT;
            z_elem = 10,
            z_min = 0,
            z_max = 1,
            radius = FT(6.371e6), # in meters
            h_elem = 10,
            n_quad_points = 4,
            staggering = CC.CommonSpaces.CellCenter(),
        )
        ones3d_level = CC.Fields.level(ones(space3d), 1)

        @test isapprox(
            Utilities.integral(ones3d_level),
            Utilities.integral(ones2d),
            rtol = 1e-5,
        )
        @test Utilities.integral(ones(space3d)) == sum(ones(space3d))
    end
end
