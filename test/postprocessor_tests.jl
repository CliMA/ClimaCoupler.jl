#=
    Unit tests for ClimaCoupler PostProcessor module
=#
import Test: @test, @testset
import ClimaCore as CC
import ClimaCoupler: PostProcessor

include("TestHelper.jl")
import .TestHelper

data(f::CC.Fields.Field) = (parent(f))
data(f::Array) = f

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")

if !Sys.iswindows() # Windows has NetCDF / HDF5 support limitations
    @testset "postprocess" begin

        cases = ((), (:regrid,), (:regrid, :zonal_mean), (:regrid, :horizontal_slice))

        results_2d = ((4, 4, 1, 96), (180, 90), (90,), (180, 90))
        results_3d = ((3, 4, 4, 1, 96), (180, 90, 3), (90, 3), (180, 90))
        expected_results = (results_2d..., results_3d...)

        ct = 0
        names = (:x,)
        for vert_nelem in (1, 3) # test (2d and 3d)
            raw_data = ones(TestHelper.create_space(Float64, nz = vert_nelem))
            for p_methods in cases
                ct += 1
                post_data = PostProcessor.postprocess(names[1], raw_data, p_methods, REGRID_DIR = REGRID_DIR)
                @test size(data(post_data.data)) == expected_results[ct]
            end

        end

    end
end
