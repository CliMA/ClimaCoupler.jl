#=
    Unit tests for ClimaCoupler Regridder module functions which require MPI

These are in a separate testing file from the other Regridder unit tests so
that MPI can be enabled for testing of these functions.
=#

using ClimaCoupler: TestHelper, Regridder
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaComms
using Dates
using Test

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")

# Set up MPI communications context
# Note that runs will hang if a context is initialized twice in the same file,
# so this context should be shared among all tests in this file.
comms_ctx = ClimaComms.SingletonCommsContext(ClimaComms.CPUDevice())
pid, nprocs = ClimaComms.init(comms_ctx)

@testset "test write_to_hdf5 and read_from_hdf5 with MPI" begin
    for FT in (Float32, Float64)
        # Set up testing directory
        mkpath(REGRID_DIR)

        hd_outfile_root = "hdf5_out_test"
        tx = Dates.DateTime(1979, 01, 01, 01, 00, 00)
        test_space = TestHelper.create_space(FT, comms_ctx = comms_ctx)
        input_field = Fields.ones(test_space)
        varname = "testdata"

        ClimaComms.barrier(comms_ctx)
        Regridder.write_to_hdf5(REGRID_DIR, hd_outfile_root, tx, input_field, varname, comms_ctx)

        output_field = Regridder.read_from_hdf5(REGRID_DIR, hd_outfile_root, tx, varname, comms_ctx)
        @test parent(input_field) == parent(output_field)
    end
end
