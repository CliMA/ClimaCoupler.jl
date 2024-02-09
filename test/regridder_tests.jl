#=
    Unit tests for ClimaCoupler Regridder module
=#

using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaComms
using Test
using NCDatasets
using Dates

using ClimaCoupler: Interfacer, Regridder, TestHelper, TimeManager, PostProcessor
import ClimaCoupler.Interfacer: get_field, name, SurfaceModelSimulation, SurfaceStub, update_field!

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")

const comms_ctx = ClimaComms.SingletonCommsContext()
const pid, nprocs = ClimaComms.init(comms_ctx)

struct TestSurfaceSimulationA <: SurfaceModelSimulation end
struct TestSurfaceSimulationB <: SurfaceModelSimulation end
struct TestSurfaceSimulationC <: SurfaceModelSimulation end
struct TestSurfaceSimulationD <: SurfaceModelSimulation end

# Initialize weights (fractions) and initial values (fields)
get_field(::TestSurfaceSimulationA, ::Val{:random}) = 1.0
get_field(::TestSurfaceSimulationB, ::Val{:random}) = 1.0
get_field(::TestSurfaceSimulationC, ::Val{:random}) = 1.0
get_field(::TestSurfaceSimulationD, ::Val{:random}) = 1.0

get_field(::TestSurfaceSimulationA, ::Val{:area_fraction}) = 0.0
get_field(::TestSurfaceSimulationB, ::Val{:area_fraction}) = 0.5
get_field(::TestSurfaceSimulationC, ::Val{:area_fraction}) = 2.0
get_field(::TestSurfaceSimulationD, ::Val{:area_fraction}) = -10.0

struct DummyStub{C} <: SurfaceModelSimulation
    cache::C
end
get_field(sim::DummyStub, ::Val{:area_fraction}) = sim.cache.area_fraction
function update_field!(sim::DummyStub, ::Val{:area_fraction}, field::Fields.Field)
    sim.cache.area_fraction .= field
end

for FT in (Float32, Float64)
    @testset "test dummmy_remap!" begin
        test_space = TestHelper.create_space(FT)
        test_field_ones = Fields.ones(test_space)
        target_field = Fields.zeros(test_space)

        Regridder.dummmy_remap!(target_field, test_field_ones)
        @test parent(target_field) == parent(test_field_ones)
    end

    @testset "test update_surface_fractions!" begin
        test_space = TestHelper.create_space(FT)
        # Construct land fraction of 0s in top half, 1s in bottom half
        land_fraction = Fields.ones(test_space)
        dims = size(parent(land_fraction))
        m = dims[1]
        n = dims[2]
        parent(land_fraction)[1:(m ÷ 2), :, :, :] .= FT(0)

        # Construct ice fraction of 0s on left, 0.5s on right
        ice_d = Fields.zeros(test_space)
        parent(ice_d)[:, (n ÷ 2 + 1):n, :, :] .= FT(0.5)

        # Construct ice fraction of 0s on left, 0.5s on right
        ocean_d = Fields.zeros(test_space)

        # Fill in only the necessary parts of the simulation
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # comms_ctx
            nothing, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # t
            Int(200), # Δt_cpl
            (; land = land_fraction, ice = Fields.zeros(test_space), ocean = Fields.zeros(test_space)), # surface_fractions
            (; ice_sim = DummyStub((; area_fraction = ice_d)), ocean_sim = SurfaceStub((; area_fraction = ocean_d))), # model_sims
            (;), # mode
            (), # diagnostics
            (;), # callbacks
            (;), # dirs
        )

        Regridder.update_surface_fractions!(cs)

        # Test that sum of fractions is 1 everywhere
        @test all(parent(cs.surface_fractions.ice .+ cs.surface_fractions.land .+ cs.surface_fractions.ocean) .== FT(1))
    end

    @testset "test combine_surfaces_from_sol!" begin
        test_space = TestHelper.create_space(FT)
        combined_field = Fields.ones(test_space)

        # Initialize weights (fractions) and initial values (fields)
        fractions = (a = 0.0, b = 0.5, c = 2.0, d = -10.0)
        fields = (a = 1.0, b = 1.0, c = 1.0, d = 1.0)

        Regridder.combine_surfaces_from_sol!(combined_field::Fields.Field, fractions::NamedTuple, fields::NamedTuple)
        @test all(parent(combined_field) .== FT(sum(fractions) * sum(fields) / length(fields)))
    end

    @testset "test combine_surfaces" begin
        test_space = TestHelper.create_space(FT)
        combined_field = Fields.ones(test_space)

        var_name = Val(:random)
        sims = (;
            a = TestSurfaceSimulationA(),
            b = TestSurfaceSimulationB(),
            c = TestSurfaceSimulationC(),
            d = TestSurfaceSimulationD(),
        )

        fractions = (
            a = get_field(sims.a, Val(:area_fraction)),
            b = get_field(sims.b, Val(:area_fraction)),
            c = get_field(sims.c, Val(:area_fraction)),
            d = get_field(sims.d, Val(:area_fraction)),
        )
        fields = (
            a = get_field(sims.a, var_name),
            b = get_field(sims.b, var_name),
            c = get_field(sims.c, var_name),
            d = get_field(sims.d, var_name),
        )

        Regridder.combine_surfaces!(combined_field, sims, var_name)
        @test all(parent(combined_field) .== FT(sum(fractions) * sum(fields) / length(fields)))
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

        @testset "test land_fraction for FT=$FT" begin
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
            land_fraction_mono =
                Regridder.land_fraction(FT, REGRID_DIR, comms_ctx, data_path, varname, test_space, mono = true)

            # Test no new extrema are introduced in monotone remapping
            nt = NCDataset(data_path) do ds
                max_val = maximum(ds[varname])
                min_val = minimum(ds[varname])
                (; max_val, min_val)
            end
            (; max_val, min_val) = nt

            @test maximum(land_fraction_mono) <= max_val
            @test minimum(land_fraction_mono) >= min_val

            # Test that monotone remapping a dataset of all ones conserves surface area
            @test sum(land_fraction_mono) - 4 * π * (R^2) < 10e-14

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
            land_fraction_halves =
                Regridder.land_fraction(FT, REGRID_DIR, comms_ctx, data_path, varname, test_space, mono = false)

            # fractioning of values below threshold should result in 0
            @test all(parent(land_fraction_halves) .== FT(0))

            # Delete testing directory and files
            rm(REGRID_DIR; recursive = true, force = true)
        end

        @testset "test hdwrite_regridfile_rll_to_cgll 3d space for FT=$FT" begin
            # Test setup
            R = FT(6371e3)
            space = TestHelper.create_space(FT, nz = 2, ne = 16, R = R)

            ispath(REGRID_DIR) ? rm(REGRID_DIR; recursive = true, force = true) : nothing
            mkpath(REGRID_DIR)

            # lat-lon dataset
            data = ones(720, 360, 2, 3) # (lon, lat, z, time)
            time = [19000101.0, 19000201.0, 19000301.0]
            lats = collect(range(-90, 90, length = 360))
            lons = collect(range(-180, 180, length = 720))
            z = [1000.0, 2000.0]
            data = reshape(sin.(lats * π / 90)[:], 1, :, 1, 1) .* data
            varname = "sinlat"

            # save the lat-lon data to a netcdf file in the required format for TempestRemap
            datafile_rll = joinpath(REGRID_DIR, "lat_lon_data.nc")
            NCDataset(datafile_rll, "c") do ds

                defDim(ds, "lat", size(lats)...)
                defDim(ds, "lon", size(lons)...)
                defDim(ds, "z", size(z)...)
                defDim(ds, "date", size(time)...)

                defVar(ds, "lon", lons, ("lon",))
                defVar(ds, "lat", lats, ("lat",))
                defVar(ds, "z", z, ("z",))
                defVar(ds, "date", time, ("date",))

                defVar(ds, varname, data, ("lon", "lat", "z", "date"))

            end

            hd_outfile_root = "data_cgll_test"
            Regridder.hdwrite_regridfile_rll_to_cgll(
                FT,
                REGRID_DIR,
                datafile_rll,
                varname,
                space,
                mono = true,
                hd_outfile_root = hd_outfile_root,
            )

            # read in data on CGLL grid from the last saved date
            date1 = TimeManager.strdate_to_datetime.(string(Int(time[end])))
            cgll_path = joinpath(REGRID_DIR, "$(hd_outfile_root)_$date1.hdf5")
            hdfreader = InputOutput.HDF5Reader(cgll_path, comms_ctx)
            T_cgll = InputOutput.read_field(hdfreader, varname)
            Base.close(hdfreader)

            # regrid back to lat-lon
            T_rll, _ = Regridder.cgll2latlonz(T_cgll)

            # check consistency across z-levels
            @test T_rll[:, :, 1] == T_rll[:, :, 2]

            # check consistency of CGLL remapped data with original data
            @test all(isapprox.(extrema(data), extrema(parent(T_cgll)), atol = 1e-2))

            # check consistency of lat-lon remapped data with original data
            @test all(isapprox.(extrema(data), extrema(T_rll), atol = 1e-3))

            # visual inspection
            # Plots.plot(T_cgll) # using ClimaCorePlots
            # Plots.contourf(Array(T_rll)[:,1])

            # Delete testing directory and files
            rm(REGRID_DIR; recursive = true, force = true)

        end
    end
end
