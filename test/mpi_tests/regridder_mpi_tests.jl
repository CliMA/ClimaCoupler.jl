#=
    Unit tests for ClimaCoupler Regridder module functions which require MPI

These are in a separate testing file so that MPI can be enabled for 
testing of these functions
=#

using ClimaCoupler, ClimaCoupler.TestHelper, ClimaCoupler.Regridder
using ClimaCore: ClimaCore, Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCommsMPI
using ClimaComms
using Dates
using Test

FT = Float64
REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")

@testset "test write_to_hdf5 and read_from_hdf5 with MPI" begin
    # Set up testing directory
    mkpath(REGRID_DIR)

    comms_ctx = ClimaCommsMPI.MPICommsContext()
    pid, nprocs = ClimaComms.init(comms_ctx)

    hd_outfile_root = "hdf5_out_test"
    tx = Dates.DateTime(1979, 01, 01, 01, 00, 00)
    test_space = create_space(FT, comms_ctx = comms_ctx)
    input_field = Fields.ones(test_space)
    varname = "testdata"

    @show comms_ctx
    ClimaComms.barrier(comms_ctx)
    write_to_hdf5(REGRID_DIR, hd_outfile_root, tx, input_field, varname, comms_ctx)

    output_field = read_from_hdf5(REGRID_DIR, hd_outfile_root, tx, varname, comms_ctx)
    @test parent(input_field) == parent(output_field)
end
