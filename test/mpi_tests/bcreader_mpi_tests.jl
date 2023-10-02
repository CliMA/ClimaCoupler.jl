#=
    Unit tests for ClimaCoupler BCReader module functions which require MPI

These are in a separate testing file from the other BCReader unit tests so
that MPI can be enabled for testing of these functions.
=#

using ClimaCoupler: Regridder, BCReader, TimeManager, Interfacer
using ClimaCore: Fields, Meshes, Domains, Topologies, Spaces
using ClimaComms
using Test
using Dates
using NCDatasets
import ArtifactWrappers as AW

# Get the path to the necessary data file - sst map
include(joinpath(@__DIR__, "..", "..", "artifacts", "artifact_funcs.jl"))
const sst_data = joinpath(sst_dataset_path(), "sst.nc")

# set up MPI communications context
device = ClimaComms.CPUSingleThreaded()
const comms_ctx = ClimaComms.context(device)
const pid, nprocs = ClimaComms.init(comms_ctx)
ClimaComms.barrier(comms_ctx)

@testset "test bcf_info_init with MPI" begin
    for FT in (Float32, Float64)
        # setup for test
        radius = FT(6731e3)
        Nq = 4
        domain = Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.DistributedTopology2D(comms_ctx, mesh, Topologies.spacefillingcurve(mesh))
        quad = Spaces.Quadratures.GLL{Nq}()
        boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)
        land_fraction_t = Fields.zeros(boundary_space_t)

        datafile_rll = sst_data
        varname = "SST"
        mono = true

        regrid_dir = "bcreader_regrid_dir"
        mkpath(regrid_dir)

        bcf_info = BCReader.bcfile_info_init(
            FT,
            regrid_dir,
            datafile_rll,
            varname,
            boundary_space_t,
            comms_ctx,
            segment_idx0 = [Int(1309)],
            land_fraction = land_fraction_t,
            mono = mono,
        )
        ClimaComms.barrier(comms_ctx)

        # test that created object exists and has correct components
        @test @isdefined(bcf_info)
        @test all(parent(bcf_info.land_fraction) .== 0)

        # construct weightfile name to test values
        hd_outfile_root = varname * "_cgll"
        outfile = hd_outfile_root * ".nc"
        outfile_root = mono ? outfile[1:(end - 3)] * "_mono" : outfile[1:(end - 3)]
        weightfile = joinpath(regrid_dir, outfile_root * "_remap_weights.nc")

        # test monotone remapping (all weights in [0, 1])
        nt = NCDataset(weightfile) do weights
            max_weight = maximum(weights["S"])
            min_weight = minimum(weights["S"])
            (; max_weight, min_weight)
        end
        (; max_weight, min_weight) = nt

        @test max_weight <= FT(1.0) || isapprox(max_weight, FT(1.0), atol = 1e-16)
        @test min_weight >= FT(0.0) || isapprox(min_weight, FT(0.0), atol = 1e-16)

        # delete testing directory and files
        ClimaComms.barrier(comms_ctx)
        ClimaComms.iamroot(comms_ctx) ? rm(regrid_dir; recursive = true, force = true) : nothing
        ClimaComms.barrier(comms_ctx)
    end
end

@testset "test update_midmonth_data! with MPI" begin
    for FT in (Float32, Float64)
        # setup for test
        date0 = date1 = DateTime(1979, 01, 01, 01, 00, 00)
        date = DateTime(1979, 01, 01, 00, 00, 00)
        tspan = (1, 90 * 86400) # Jan-Mar
        Δt = 1 * 3600

        radius = FT(6731e3)
        Nq = 4
        domain = Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.DistributedTopology2D(comms_ctx, mesh, Topologies.spacefillingcurve(mesh))
        quad = Spaces.Quadratures.GLL{Nq}()
        boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)

        land_fraction_t = Fields.zeros(boundary_space_t)
        dummy_data = (; test_data = zeros(axes(land_fraction_t)))

        datafile_rll = sst_data
        varname = "SST"

        regrid_dir = "bcreader_regrid_dir"
        mkpath(regrid_dir)

        ClimaComms.barrier(comms_ctx)

        bcf_info = BCReader.bcfile_info_init(
            FT,
            regrid_dir,
            datafile_rll,
            varname,
            boundary_space_t,
            comms_ctx,
            segment_idx0 = [Int(1309)],
            interpolate_daily = false,
            land_fraction = land_fraction_t,
        )

        dates = (; date = [date], date0 = [date0], date1 = [date1])
        SST_all = []
        updating_dates = []

        cs_t = Interfacer.CoupledSimulation{FT}(
            comms_ctx, # comms_ctx
            dates, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            tspan, # tspan
            Int(0), # t
            Δt, # Δt_cpl
            (;), # surface_fractions
            (;), # model_sims
            (;), # mode
            (), # diagnostics
        )

        ClimaComms.barrier(comms_ctx)

        # step in time
        walltime = @elapsed for t in ((tspan[1] + Δt):Δt:tspan[end])
            cs_t.dates.date[1] = TimeManager.current_date(cs_t, t) # if not global, `date`` is not updated. Check that this still runs when distributed.

            model_date = cs_t.dates.date[1]
            callback_date = BCReader.next_date_in_file(bcf_info)

            # TODO investigate if macro would be faster here
            if (model_date >= callback_date)
                BCReader.update_midmonth_data!(model_date, bcf_info)
                push!(SST_all, deepcopy(bcf_info.monthly_fields[1]))
                push!(updating_dates, deepcopy(model_date))
            end

        end

        ClimaComms.barrier(comms_ctx)
        # test if the SST field was modified
        @test SST_all[end] !== SST_all[end - 1]
        # check that the final file date is as expected
        @test Date(updating_dates[end]) == Date(1979, 03, 16)

        # test warning/error cases
        current_fields = Fields.zeros(FT, boundary_space_t), Fields.zeros(FT, boundary_space_t)

        # use this function to reset values between test cases
        function reset_bcf_info(bcf_info)
            bcf_info.monthly_fields[1] .= current_fields[1]
            bcf_info.monthly_fields[2] .= current_fields[2]
            bcf_info.segment_length[1] = Int(1)
        end

        hd_outfile_root = varname * "_cgll"

        #  case 1: segment_idx == segment_idx0, date < all_dates[segment_idx]
        bcf_info.segment_idx[1] = bcf_info.segment_idx0[1]
        date = DateTime(bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1))
        BCReader.update_midmonth_data!(date, bcf_info)

        ClimaComms.barrier(comms_ctx)
        @test bcf_info.monthly_fields[1] == bcf_info.scaling_function(
            Regridder.read_from_hdf5(
                regrid_dir,
                hd_outfile_root,
                bcf_info.all_dates[Int(bcf_info.segment_idx0[1])],
                varname,
                comms_ctx,
            ),
            bcf_info,
        )
        @test bcf_info.monthly_fields[2] == bcf_info.monthly_fields[1]
        @test bcf_info.segment_length[1] == Int(0)

        #  case 2: date > all_dates[end - 1]
        reset_bcf_info(bcf_info)
        date = DateTime(bcf_info.all_dates[end - 1] + Dates.Day(1))
        BCReader.update_midmonth_data!(date, bcf_info)

        ClimaComms.barrier(comms_ctx)

        @test bcf_info.monthly_fields[1] == bcf_info.scaling_function(
            Regridder.read_from_hdf5(
                regrid_dir,
                hd_outfile_root,
                bcf_info.all_dates[Int(length(bcf_info.all_dates))],
                varname,
                comms_ctx,
            ),
            bcf_info,
        )
        @test bcf_info.monthly_fields[2] == bcf_info.monthly_fields[1]
        @test bcf_info.segment_length[1] == Int(0)

        #  case 3: date - all_dates[segment_idx + 1] > 2
        reset_bcf_info(bcf_info)
        date = DateTime(bcf_info.all_dates[bcf_info.segment_idx[1] + 1] + Dates.Day(3))
        BCReader.update_midmonth_data!(date, bcf_info)

        nearest_idx = argmin(
            abs.(
                parse(FT, TimeManager.datetime_to_strdate(date)) .-
                parse.(FT, TimeManager.datetime_to_strdate.(bcf_info.all_dates[:]))
            ),
        )

        ClimaComms.barrier(comms_ctx)
        @test bcf_info.segment_idx[1] == bcf_info.segment_idx0[1] == nearest_idx

        #  case 4: everything else
        reset_bcf_info(bcf_info)
        bcf_info.segment_idx[1] = bcf_info.segment_idx0[1] + Int(1)
        date = bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1)

        ClimaComms.barrier(comms_ctx)
        @test_throws ErrorException BCReader.update_midmonth_data!(date, bcf_info)

        ClimaComms.barrier(comms_ctx)
        ClimaComms.iamroot(comms_ctx) ? rm(regrid_dir; recursive = true, force = true) : nothing
        ClimaComms.barrier(comms_ctx)
    end
end
