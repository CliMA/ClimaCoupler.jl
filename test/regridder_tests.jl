#=
    Unit tests for ClimaCoupler Regridder module
=#

using ClimaCoupler: Utilities, Regridder, TestHelper
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaComms
using Test
using NCDatasets
using Dates

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")

const comms_ctx = ClimaComms.SingletonCommsContext()
const pid, nprocs = ClimaComms.init(comms_ctx)

for FT in (Float32, Float64)
    @testset "test dummmy_remap!" begin
        test_space = TestHelper.create_space(FT)
        test_field_ones = Fields.ones(test_space)
        target_field = Fields.zeros(test_space)

        Regridder.dummmy_remap!(target_field, test_field_ones)
        @test parent(target_field) == parent(test_field_ones)
    end

    @testset "test update_masks!" begin
        test_space = TestHelper.create_space(FT)
        # Construct land mask of 0s in top half, 1s in bottom half
        land_mask = Fields.ones(test_space)
        dims = size(parent(land_mask))
        m = dims[1]
        n = dims[2]
        parent(land_mask)[1:(m ÷ 2), :, :, :] .= FT(0)

        # Construct ice mask of 0s on left, 0.5s on right
        ice_d = Fields.zeros(test_space)
        parent(ice_d)[:, (n ÷ 2 + 1):n, :, :] .= FT(0.5)

        # Fill in only the necessary parts of the simulation
        cs = Utilities.CoupledSimulation{FT}(
            nothing, # comms_ctx
            nothing, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # t
            Int(200), # Δt_cpl
            (; land = land_mask, ice = Fields.zeros(test_space), ocean = Fields.zeros(test_space)), # surface_masks
            (; ice_sim = (; integrator = (; p = (; ice_fraction = ice_d)))), # model_sims
            (;), # mode
            (), # diagnostics
        )

        Regridder.update_masks!(cs)

        # Test that sum of masks is 1 everywhere
        @test all(parent(cs.surface_masks.ice .+ cs.surface_masks.land .+ cs.surface_masks.ocean) .== FT(1))
    end

    @testset "test combine_surfaces!" begin
        test_space = TestHelper.create_space(FT)
        combined_field = Fields.ones(test_space)

        # Initialize weights (masks) and initial values (fields)
        masks = (a = 0.0, b = 0.5, c = 2.0, d = -10.0)
        fields = (a = 1.0, b = 1.0, c = 1.0, d = 1.0)

        Regridder.combine_surfaces!(combined_field::Fields.Field, masks::NamedTuple, fields::NamedTuple)
        @test all(parent(combined_field) .== FT(sum(masks) * sum(fields) / length(fields)))
    end

    # Add tests which use TempestRemap here -
    # TempestRemap is not built on Windows because of NetCDF support limitations
    if !Sys.iswindows()
        @testset "test write_to_hdf5 and read_from_hdf5" begin
            # Set up testing directory
            ispath(REGRID_DIR) ? rm(REGRID_DIR; recursive = true, force = true) : nothing
            mkpath(REGRID_DIR)

            hd_outfile_root = "hdf5_out_test"
            tx = Dates.DateTime(1979, 01, 01, 01, 00, 00)
            test_space = TestHelper.create_space(FT)
            input_field = Fields.ones(test_space)
            varname = "testdata"

            Regridder.write_to_hdf5(REGRID_DIR, hd_outfile_root, tx, input_field, varname, comms_ctx)

            output_field = Regridder.read_from_hdf5(REGRID_DIR, hd_outfile_root, tx, varname, comms_ctx)
            @test parent(input_field) == parent(output_field)

            # Delete testing directory and files
            rm(REGRID_DIR; recursive = true, force = true)
        end

        @testset "test remap_field_cgll_to_rll for FT=$FT" begin
            # Set up testing directory
            remap_tmpdir = joinpath(REGRID_DIR, "cgll_to_rll")
            ispath(remap_tmpdir) ? rm(remap_tmpdir; recursive = true, force = true) : nothing
            mkpath(remap_tmpdir)
            name = "testdata"
            datafile_rll = remap_tmpdir * "/" * name * "_rll.nc"

            test_space = TestHelper.create_space(FT)
            field = Fields.ones(test_space)

            Regridder.remap_field_cgll_to_rll(name, field, remap_tmpdir, datafile_rll)

            # Test no new extrema are introduced in monotone remapping
            nt = NCDataset(datafile_rll) do ds
                max_remapped = maximum(ds[name])
                min_remapped = minimum(ds[name])
                (; max_remapped, min_remapped)
            end
            (; max_remapped, min_remapped) = nt

            @test max_remapped <= maximum(field)
            @test min_remapped >= minimum(field)

            # Delete testing directory and files
            rm(REGRID_DIR; recursive = true, force = true)
        end

        @testset "test land_sea_mask for FT=$FT" begin
            # Test setup
            R = FT(6371e3)
            test_space = TestHelper.create_space(FT, R = R)
            ispath(REGRID_DIR) ? rm(REGRID_DIR; recursive = true, force = true) : nothing
            mkpath(REGRID_DIR)

            # Initialize dataset of all ones
            data_path = joinpath(REGRID_DIR, "ls_mask_data.nc")
            varname = "test_data"
            TestHelper.gen_ncdata(FT, data_path, varname, FT(1))

            # Test monotone masking
            land_mask_mono =
                Regridder.land_sea_mask(FT, REGRID_DIR, comms_ctx, data_path, varname, test_space, mono = true)

            # Test no new extrema are introduced in monotone remapping
            nt = NCDataset(data_path) do ds
                max_val = maximum(ds[varname])
                min_val = minimum(ds[varname])
                (; max_val, min_val)
            end
            (; max_val, min_val) = nt

            @test maximum(land_mask_mono) <= max_val
            @test minimum(land_mask_mono) >= min_val

            # Test that monotone remapping a dataset of all ones conserves surface area
            @test sum(land_mask_mono) - 4 * π * (R^2) < 10e-14

            # Delete testing directory and files
            rm(REGRID_DIR; recursive = true, force = true)

            # Set up testing directory
            ispath(REGRID_DIR) ? rm(REGRID_DIR; recursive = true, force = true) : nothing
            mkpath(REGRID_DIR)

            # Initialize dataset of all 0.5s
            data_path = joinpath(REGRID_DIR, "ls_mask_data.nc")
            varname = "test_data_halves"
            TestHelper.gen_ncdata(FT, data_path, varname, FT(0.5))

            # Test non-monotone masking
            land_mask_halves =
                Regridder.land_sea_mask(FT, REGRID_DIR, comms_ctx, data_path, varname, test_space, mono = false)

            # Masking of values below threshold should result in 0
            @test all(parent(land_mask_halves) .== FT(0))

            # Delete testing directory and files
            rm(REGRID_DIR; recursive = true, force = true)
        end
    end
end
