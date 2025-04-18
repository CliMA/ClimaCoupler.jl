#=
    Unit tests for ClimaCoupler Utilities module
=#
import Test: @testset, @test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCoupler: Utilities
import ClimaCore as CC

# Initialize MPI context, in case
ClimaComms.init(ClimaComms.context())

for FT in (Float32, Float64)
    @testset "test swap_space!" begin
        space1 = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4)
        space2 = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4)

        field1 = ones(space1)
        field2 = ones(space2)

        field2 = Utilities.swap_space!(space2, field1)

        @test parent(field1) == parent(field2)
        @test axes(field2) == space2
    end

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

    @testset "test binary_mask" begin
        space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4)
        @test all(parent(Utilities.binary_mask.(zeros(space))) .== 0)
        @test all(parent(Utilities.binary_mask.(ones(space))) .== 1)
        @test all(parent(Utilities.binary_mask.(fill(FT(0.5), space), FT(0.6))) .== 0)
        @test all(parent(Utilities.binary_mask.(fill(FT(0.5), space), FT(0.4))) .== 1)
    end
end
