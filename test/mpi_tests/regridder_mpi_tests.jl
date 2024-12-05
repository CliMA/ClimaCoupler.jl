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
import ClimaUtilities.Regridders
import ClimaCoreTempestRemap

include(joinpath("..", "TestHelper.jl"))
import .TestHelper

# Set up MPI communications context
# Note that runs will hang if a context is initialized twice in the same file,
# so this context should be shared among all tests in this file.
device = ClimaComms.CPUSingleThreaded()
const comms_ctx = ClimaComms.context(device)
pid, nprocs = ClimaComms.init(comms_ctx)

@testset "test read_from_hdf5 with MPI" begin
    for FT in (Float32, Float64)
        # Create temporary regrid directory on root process and broadcast
        regrid_dir = nothing
        if ClimaComms.iamroot(comms_ctx)
            regrid_dir = mktempdir(pwd(), prefix = "regrid_tmp_")
        end
        regrid_dir = ClimaComms.bcast(comms_ctx, regrid_dir)

        test_space = TestHelper.create_space(FT, comms_ctx = comms_ctx)
        varname = "testdata"
        hd_outfile_root = varname
        data_path = joinpath(regrid_dir, "test_data")
        TestHelper.gen_ncdata_time(FT, data_path, varname, FT(1))

        ClimaComms.barrier(comms_ctx)
        Regridders.TempestRegridder(test_space, varname, data_path; regrid_dir, mono = true)
        ClimaComms.barrier(comms_ctx)
        output_field_ones =
            Regridder.read_from_hdf5(regrid_dir, hd_outfile_root, Dates.DateTime(2021), varname, comms_ctx)
        output_field_zeros =
            Regridder.read_from_hdf5(regrid_dir, hd_outfile_root, Dates.DateTime(2020), varname, comms_ctx)
        @test parent(CC.Fields.ones(test_space)) == parent(output_field_ones)
        @test parent(CC.Fields.zeros(test_space)) == parent(output_field_zeros)
    end
end
