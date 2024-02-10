#=
    Unit tests for ClimaCoupler BCReader module
=#

using ClimaCoupler: Regridder, BCReader, TimeManager, Interfacer
using ClimaCore: Fields, Meshes, Domains, Topologies, Spaces
using ClimaComms
using Test
using Dates
using NCDatasets
import ArtifactWrappers as AW

# get the paths to the necessary data files - sst map, land sea mask
include(joinpath(@__DIR__, "..", "artifacts", "artifact_funcs.jl"))
sst_data = joinpath(sst_dataset_path(), "sst.nc")
mask_data = joinpath(mask_dataset_path(), "seamask.nc")

const comms_ctx = ClimaComms.SingletonCommsContext()
const pid, nprocs = ClimaComms.init(comms_ctx)

for FT in (Float32, Float64)
    @testset "test next_date_in_file for FT=$FT" begin
        dummy_dates = Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 10))
        date0 = dummy_dates[1]
        segment_idx0 = [
            argmin(
                abs.(
                    parse(FT, TimeManager.datetime_to_strdate(date0)) .-
                    parse.(FT, TimeManager.datetime_to_strdate.(dummy_dates[:]))
                ),
            ),
        ]

        bcf_info = BCReader.BCFileInfo{FT}(
            "",                                 # bcfile_dir
            comms_ctx,                          # comms_ctx
            "",                                 # hd_outfile_root
            "",                                 # varname
            dummy_dates,                        # all_dates
            nothing,                            # monthly_fields
            nothing,                            # scaling_function
            nothing,                            # land_fraction
            deepcopy(segment_idx0),             # segment_idx
            segment_idx0,                       # segment_idx0
            Int[],                              # segment_length
            false,                              # interpolate_daily
        )

        idx = segment_idx0[1]
        current_date = date0
        next_date = BCReader.next_date_in_file(bcf_info)
        @test current_date == dummy_dates[idx]

        for i in 1:(length(dummy_dates) - 2)
            current_date = next_date
            bcf_info.segment_idx[1] += Int(1)
            next_date = BCReader.next_date_in_file(bcf_info)
            idx = segment_idx0[1] + i
            @test next_date == dummy_dates[idx + 1]
        end
    end

    @testset "test interpol for FT=$FT" begin
        # Setup
        t1 = FT(0)
        t2 = FT(10)

        f1 = FT(-1)
        f2 = FT(1)

        # Case 1: t1 < t < t2 --> should return average (0)
        t = FT(5)
        @test BCReader.interpol(f1, f2, t - t1, t2 - t1) == 0

        # Case 2: t1 = t < t2 --> should return val at t1 (f1)
        t = FT(0)
        @test BCReader.interpol(f1, f2, t - t1, t2 - t1) == f1

        # Case 3: t1 < t = t2 --> should return val at t2 (f2)
        t = FT(10)
        @test BCReader.interpol(f1, f2, t - t1, t2 - t1) == f2
    end

    @testset "test interpolate_midmonth_to_daily for FT=$FT" begin
        # test interpolate_midmonth_to_daily with interpolation
        interpolate_daily = true
        dummy_dates = Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 100))
        segment_idx0 = [Int(1)]

        # these values give an `interp_fraction` of 0.5 in `interpol` for ease of testing
        date0 = dummy_dates[Int(segment_idx0[1] + 1)]
        segment_length = [Int(2) * ((date0 - dummy_dates[Int(segment_idx0[1])]).value)]

        radius = FT(6731e3)
        Nq = 4
        domain = Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.Topology2D(comms_ctx, mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)
        monthly_fields = (zeros(boundary_space_t), ones(boundary_space_t))

        bcf_info_interp = BCReader.BCFileInfo{FT}(
            "",                                 # bcfile_dir
            comms_ctx,                          # comms_ctx
            "",                                 # hd_outfile_root
            "",                                 # varname
            dummy_dates,                        # all_dates
            monthly_fields,                     # monthly_fields
            nothing,                            # scaling_function
            nothing,                            # land_fraction
            deepcopy(segment_idx0),             # segment_idx
            segment_idx0,                       # segment_idx0
            segment_length,                     # segment_length
            interpolate_daily,                  # interpolate_daily
        )
        @test BCReader.interpolate_midmonth_to_daily(date0, bcf_info_interp) == ones(boundary_space_t) .* FT(0.5)

        # test interpolate_midmonth_to_daily without interpolation
        interpolate_daily = false

        bcf_info_no_interp = BCReader.BCFileInfo{FT}(
            "",                                 # bcfile_dir
            comms_ctx,                          # comms_ctx
            "",                                 # hd_outfile_root
            "",                                 # varname
            dummy_dates,                        # all_dates
            monthly_fields,                     # monthly_fields
            nothing,                            # scaling_function
            nothing,                            # land_fraction
            deepcopy(segment_idx0),             # segment_idx
            segment_idx0,                       # segment_idx0
            segment_length,                     # segment_length
            interpolate_daily,                  # interpolate_daily
        )
        @test BCReader.interpolate_midmonth_to_daily(date0, bcf_info_no_interp) == monthly_fields[1]
    end

    # Add tests which use TempestRemap here -
    # TempestRemap is not built on Windows because of NetCDF support limitations
    # `bcf_info_init` uses TR via a call to `hdwrite_regridfile_rll_to_cgll`
    if !Sys.iswindows()
        @testset "test update_midmonth_data! for FT=$FT" begin
            # setup for test
            date0 = date1 = DateTime(1979, 01, 01, 01, 00, 00)
            date = DateTime(1979, 01, 01, 00, 00, 00)
            tspan = (Int(1), Int(90 * 86400)) # Jan-Mar
            Δt = Int(1 * 3600)

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
            isdir(regrid_dir) ? nothing : mkpath(regrid_dir)

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
            )

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
                bcf_info.segment_idx[1] = bcf_info.segment_idx0[1]
            end

            hd_outfile_root = varname * "_cgll"

            #  case 1: date < all_dates[segment_idx] (init)
            bcf_info.segment_idx[1] = bcf_info.segment_idx0[1]
            date = DateTime(bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1))
            BCReader.update_midmonth_data!(date, bcf_info)

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
                date = DateTime(bcf_info.all_dates[bcf_info.segment_idx0[1]]) + extra
                BCReader.update_midmonth_data!(date, bcf_info)

                end_field_c2 = deepcopy(bcf_info.monthly_fields[2])
                segment_length_c2 = deepcopy(bcf_info.segment_length[1])
                current_index_c2 = deepcopy(bcf_info.segment_idx[1])

                # modified field
                @test end_field_c2 !== bcf_info.monthly_fields[1]
                # updated segment length
                @test segment_length_c2[1] !== Int(0)
                # updated reference segment index
                @test current_index_c2 == bcf_info.segment_idx0[1] + 1

                # case 2: (date - all_dates[Int(segment_idx0) + 1]) >= 0 (init)
                # do not reset segment_idx0. It's current value ensures that we get the same result as case 3
                reset_bcf_info(bcf_info)

                date = DateTime(bcf_info.all_dates[bcf_info.segment_idx0[1] + 1]) + extra
                BCReader.update_midmonth_data!(date, bcf_info)

                nearest_idx = argmin(
                    abs.(
                        parse(FT, TimeManager.datetime_to_strdate(date)) .-
                        parse.(FT, TimeManager.datetime_to_strdate.(bcf_info.all_dates[:]))
                    ),
                )

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

                date = DateTime(bcf_info.all_dates[bcf_info.segment_idx0[1]]) + extra
                BCReader.update_midmonth_data!(date, bcf_info)

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

                date = DateTime(bcf_info.all_dates[bcf_info.segment_idx0[1]] + extra)
                BCReader.update_midmonth_data!(date, bcf_info)

                @test bcf_info.segment_idx[1] == bcf_info.segment_idx0[1] + 1
            end

            #  case 6: everything else
            reset_bcf_info(bcf_info)
            bcf_info.segment_idx[1] = bcf_info.segment_idx0[1] + Int(1)
            date = bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1)

            @test_throws ErrorException BCReader.update_midmonth_data!(date, bcf_info)

            rm(regrid_dir; recursive = true, force = true)
        end

        @testset "test bcf_info_init for FT=$FT" begin
            # setup for test
            radius = FT(6731e3)
            Nq = 4
            domain = Domains.SphereDomain(radius)
            mesh = Meshes.EquiangularCubedSphere(domain, 4)
            topology = Topologies.Topology2D(comms_ctx, mesh)
            quad = Spaces.Quadratures.GLL{Nq}()
            boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)
            land_fraction_t = Fields.zeros(boundary_space_t)

            datafile_rll = mask_data
            varname = "LSMASK"
            mono = true

            regrid_dir = "bcreader_regrid_dir"
            isdir(regrid_dir) ? nothing : mkpath(regrid_dir)

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
            rm(regrid_dir; recursive = true, force = true)
        end
    end
end
