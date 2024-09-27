#=
    Unit tests for ClimaCoupler Regridder module functions which require MPI

These are in a separate testing file from the other Regridder unit tests so
that MPI can be enabled for testing of these functions.
=#
import Test: @test, @testset
import Dates
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaCoupler
import ClimaCoupler: Regridder

include(joinpath("..", "TestHelper.jl"))
import .TestHelper

# Set up MPI communications context
# Note that runs will hang if a context is initialized twice in the same file,
# so this context should be shared among all tests in this file.
device = ClimaComms.CPUSingleThreaded()
const comms_ctx = ClimaComms.context(device)
pid, nprocs = ClimaComms.init(comms_ctx)

@testset "test write_to_hdf5 and read_from_hdf5 with MPI" begin
    for FT in (Float32, Float64)
        # Create temporary regrid directory on root process and broadcast
        regrid_dir = nothing
        if ClimaComms.iamroot(comms_ctx)
            regrid_dir = mktempdir(pwd(), prefix = "regrid_tmp_")
        end
        regrid_dir = ClimaComms.bcast(comms_ctx, regrid_dir)

        hd_outfile_root = "hdf5_out_test"
        tx = Dates.DateTime(1979, 01, 01, 01, 00, 00)
        test_space = TestHelper.create_space(FT, comms_ctx = comms_ctx)
        input_field = CC.Fields.ones(test_space)
        varname = "testdata"

        ClimaComms.barrier(comms_ctx)
        Regridder.write_to_hdf5(regrid_dir, hd_outfile_root, tx, input_field, varname, comms_ctx)

        ClimaComms.barrier(comms_ctx)
        output_field = Regridder.read_from_hdf5(regrid_dir, hd_outfile_root, tx, varname, comms_ctx)
        @test parent(input_field) == parent(output_field)
    end
end
