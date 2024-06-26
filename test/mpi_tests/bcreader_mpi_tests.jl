#=
    Unit tests for ClimaCoupler BCReader module functions which require MPI

These are in a separate testing file from the other BCReader unit tests so
that MPI can be enabled for testing of these functions.
=#
import Test: @test, @testset, @test_throws
import Dates
import NCDatasets
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaCoupler
import ClimaCoupler: Regridder, BCReader, TimeManager, Interfacer

# Get the path to the necessary data file - sst map
pkg_dir = pkgdir(ClimaCoupler)
include(joinpath(pkg_dir, "artifacts", "artifact_funcs.jl"))
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
        domain = CC.Domains.SphereDomain(radius)
        mesh = CC.Meshes.EquiangularCubedSphere(domain, 4)
        topology = CC.Topologies.DistributedTopology2D(comms_ctx, mesh, CC.Topologies.spacefillingcurve(mesh))
        quad = CC.Spaces.Quadratures.GLL{Nq}()
        boundary_space_t = CC.Spaces.SpectralElementSpace2D(topology, quad)
        land_fraction_t = CC.Fields.zeros(boundary_space_t)

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
        nt = NCDatasets.NCDataset(weightfile) do weights
            max_weight = maximum(weights["S"])
            min_weight = minimum(weights["S"])
            (; max_weight, min_weight)
        end
        (; max_weight, min_weight) = nt

        @test max_weight <= FT(1.0) || isapprox(max_weight, FT(1.0), atol = 1e-16)
        @test min_weight >= FT(0.0) || isapprox(min_weight, FT(0.0), atol = 1e-16)

        # delete testing directory and files
        ClimaComms.barrier(comms_ctx)
        ClimaComms.iamroot(comms_ctx) && rm(regrid_dir; recursive = true, force = true)
        ClimaComms.barrier(comms_ctx)
    end
end

@testset "test update_midmonth_data! with MPI" begin
    for FT in (Float32, Float64)
        # setup for test
        date0 = date1 = Dates.DateTime(1979, 01, 01, 01, 00, 00)
        date = Dates.DateTime(1979, 01, 01, 00, 00, 00)
        tspan = (1, 90 * 86400) # Jan-Mar
        Δt = 1 * 3600

        radius = FT(6731e3)
        Nq = 4
        domain = CC.Domains.SphereDomain(radius)
        mesh = CC.Meshes.EquiangularCubedSphere(domain, 4)
        topology = CC.Topologies.DistributedTopology2D(comms_ctx, mesh, CC.Topologies.spacefillingcurve(mesh))
        quad = CC.Spaces.Quadratures.GLL{Nq}()
        boundary_space_t = CC.Spaces.SpectralElementSpace2D(topology, quad)

        land_fraction_t = CC.Fields.zeros(boundary_space_t)
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
            (;), # callbacks
            (;), # dirs
            nothing, # turbulent_fluxes
            nothing, # thermo_params
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
        @test Dates.Date(updating_dates[end]) == Dates.Date(1979, 03, 16)

        # test warning/error cases
        current_fields = CC.Fields.zeros(FT, boundary_space_t), CC.Fields.zeros(FT, boundary_space_t)

        # use this function to reset values between test cases
        function reset_bcf_info(bcf_info)
            bcf_info.monthly_fields[1] .= current_fields[1]
            bcf_info.monthly_fields[2] .= current_fields[2]
            bcf_info.segment_length[1] = Int(1)
            bcf_info.segment_idx[1] = bcf_info.segment_idx0[1]
        end

        hd_outfile_root = varname * "_cgll"

        #  case 1: date < all_dates[segment_idx] (init)
        bcf_info.segment_idx[1] = bcf_info.segment_idx0[1]
        date = Dates.DateTime(bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1))
        BCReader.update_midmonth_data!(date, bcf_info)

        ClimaComms.barrier(comms_ctx)
        # unmodified field
        @test bcf_info.monthly_fields[2] == bcf_info.monthly_fields[1]
        # zero segment length
        @test bcf_info.segment_length[1] == Int(0)
        # segment index is reset
        @test bcf_info.segment_idx0[1] == bcf_info.segment_idx[1] - 1

        # cases 2 and 3
        extra_days = [Dates.Day(0), Dates.Day(3)]
        for extra in extra_days
            # case 3: (date - all_dates[Int(segment_idx0)]) >= 0 (init)
            reset_bcf_info(bcf_info)
            date = Dates.DateTime(bcf_info.all_dates[bcf_info.segment_idx0[1]]) + extra
            BCReader.update_midmonth_data!(date, bcf_info)

            end_field_c2 = deepcopy(bcf_info.monthly_fields[2])
            segment_length_c2 = deepcopy(bcf_info.segment_length[1])
            current_index_c2 = deepcopy(bcf_info.segment_idx[1])

            ClimaComms.barrier(comms_ctx)
            # modified field
            @test end_field_c2 !== bcf_info.monthly_fields[1]
            # updated segment length
            @test segment_length_c2[1] !== Int(0)
            # updated reference segment index
            @test current_index_c2 == bcf_info.segment_idx0[1] + 1

            # case 2: (date - all_dates[Int(segment_idx0) + 1]) >= 0 (init)
            # do not reset segment_idx0. It's current value ensures that we get the same result as case 3
            reset_bcf_info(bcf_info)

            date = Dates.DateTime(bcf_info.all_dates[bcf_info.segment_idx0[1] + 1]) + extra
            BCReader.update_midmonth_data!(date, bcf_info)

            nearest_idx = argmin(
                abs.(
                    parse(FT, TimeManager.datetime_to_strdate(date)) .-
                    parse.(FT, TimeManager.datetime_to_strdate.(bcf_info.all_dates[:]))
                ),
            )

            ClimaComms.barrier(comms_ctx)
            @test bcf_info.segment_idx[1] == bcf_info.segment_idx0[1] + 1 == nearest_idx

            # compare to case 3 (expecting the same result - this defaults to it):
            @test bcf_info.monthly_fields[1] == bcf_info.scaling_function(
                Regridder.read_from_hdf5(
                    regrid_dir,
                    hd_outfile_root,
                    bcf_info.all_dates[Int(bcf_info.segment_idx[1])],
                    varname,
                    comms_ctx,
                ),
                bcf_info,
            )

            # check case 2 defaults to case 3
            @test end_field_c2 !== bcf_info.monthly_fields[1]
            @test segment_length_c2[1] !== Int(0)
            @test current_index_c2 == bcf_info.segment_idx0[1] + 1

        end

        #  case 4: date > all_dates[end]
        for extra in extra_days
            bcf_info.segment_idx0[1] = length(bcf_info.all_dates)
            reset_bcf_info(bcf_info)

            date = Dates.DateTime(bcf_info.all_dates[bcf_info.segment_idx0[1]]) + extra
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
        end

        #  case 5: Dates.days(date - all_dates[segment_idx]) >= 0

        extra = Dates.Day(3)
        for extra in extra_days
            bcf_info.segment_idx0[1] = 2
            reset_bcf_info(bcf_info)

            date = Dates.DateTime(bcf_info.all_dates[bcf_info.segment_idx0[1]] + extra)
            BCReader.update_midmonth_data!(date, bcf_info)

            ClimaComms.barrier(comms_ctx)
            @test bcf_info.segment_idx[1] == bcf_info.segment_idx0[1] + 1
        end

        #  case 6: everything else
        reset_bcf_info(bcf_info)
        bcf_info.segment_idx[1] = bcf_info.segment_idx0[1] + Int(1)
        date = bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1)

        ClimaComms.barrier(comms_ctx)
        @test_throws ErrorException BCReader.update_midmonth_data!(date, bcf_info)

        ClimaComms.barrier(comms_ctx)
        ClimaComms.iamroot(comms_ctx) && rm(regrid_dir; recursive = true, force = true)
        ClimaComms.barrier(comms_ctx)
    end
end
