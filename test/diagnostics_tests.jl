#=
    Unit tests for ClimaCoupler Diagnostics module
=#
import Test: @test, @testset
import Dates
import ClimaComms
import ClimaCore as CC
import ClimaCoupler: ConservationChecker, Diagnostics, Interfacer, TestHelper, TimeManager

Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:x}) = 1

for FT in (Float32, Float64)
    @testset "init_diagnostics for FT=$FT" begin
        names = (:x, :y)
        space = TestHelper.create_space(FT)
        dg = Diagnostics.init_diagnostics(names, space)
        @test typeof(dg) == Diagnostics.DiagnosticsGroup{TimeManager.EveryTimestep, NamedTuple{(), Tuple{}}}
    end

    @testset "accumulate_diagnostics!, collect_diags, iterate_operations, operation{accumulation{TimeMean, Nothing}}, get_var for FT=$FT" begin
        cases = (nothing, Diagnostics.TimeMean([Int(0)]))
        expected_results = (FT(2), FT(3))
        for (c_i, case) in enumerate(cases)
            names = (:x,)
            space = TestHelper.create_space(FT)
            dg_2d = Diagnostics.init_diagnostics(
                names,
                space,
                save = TimeManager.EveryTimestep(),
                operations = (; accumulate = case),
            )
            dg_2d.field_vector .= FT(2)
            cs = Interfacer.CoupledSimulation{FT}(
                nothing, # comms_ctx
                nothing, # dates
                nothing, # boundary_space
                nothing, # fields
                nothing, # parsed_args
                nothing, # conservation_checks
                (Int(0), Int(1000)), # tspan
                Int(100), # t
                Int(100), # Δt_cpl
                (;), # surface_masks
                (;), # model_sims
                (;), # mode
                (dg_2d,),
                (;), # callbacks
                (;), # dirs
                nothing, # turbulent_fluxes
                nothing, # thermo_params
            )
            Diagnostics.accumulate_diagnostics!(cs)
            @test cs.diagnostics[1].field_vector[1] == expected_results[c_i]

            @test isnothing(Diagnostics.get_var(cs, Val(:z)))
        end
    end

    if !Sys.iswindows() # Windows has NetCDF / HDF5 support limitations
        @testset "save_diagnostics" begin
            test_dir = "diag_test_dir"
            names = (:x,)
            space = TestHelper.create_space(FT)
            dg_2d = Diagnostics.init_diagnostics(
                names,
                space,
                save = TimeManager.EveryTimestep(),
                operations = (; accumulate = Diagnostics.TimeMean([Int(0)])),
                output_dir = test_dir,
            ) # or use accumulate = nothing for snapshop save
            cs = Interfacer.CoupledSimulation{FT}(
                ClimaComms.SingletonCommsContext(), # comms_ctx
                (date = [Dates.DateTime(0, 2)], date1 = [Dates.DateTime(0, 1)]), # dates
                nothing, # boundary_space
                nothing, # fields
                nothing, # parsed_args
                nothing, # conservation_checks
                (Int(0), Int(1000)),# tspan
                Int(100), # t
                Int(100), # Δt_cpl
                (;), # surface_masks
                (;), # model_sims
                (;), # mode
                (dg_2d,), # diagnostics
                (;), # callbacks
                (;), # dirs
                nothing, # turbulent_fluxes
                nothing, # thermo_params
            )
            Diagnostics.save_diagnostics(cs, cs.diagnostics[1])
            file = filter(x -> endswith(x, ".hdf5"), readdir(test_dir))
            @test !isempty(file)
            rm(test_dir; recursive = true, force = true)

        end
    end

    @testset "save_time_format for FT=$FT" begin
        date = Dates.DateTime(1970, 2, 1, 0, 1)
        unix = Diagnostics.save_time_format(date, TimeManager.Monthly())
        @test unix == 0
    end

    @testset "pre_save{TimeMean, Nothing}, post_save for FT=$FT" begin
        cases = (nothing, Diagnostics.TimeMean([Int(0)]))
        expected_results = ((FT(3), FT(1), FT(1)), (FT(4), FT(2.5), FT(0)))

        for (c_i, case) in enumerate(cases)
            names = (:x,)
            space = TestHelper.create_space(FT)
            dg_2d = Diagnostics.init_diagnostics(
                names,
                space,
                save = TimeManager.EveryTimestep(),
                operations = (; accumulate = case),
            )
            dg_2d.field_vector .= FT(3)
            cs = Interfacer.CoupledSimulation{FT}(
                nothing, # comms_ctx
                nothing, # dates
                nothing, # boundary_space
                nothing, # fields
                nothing, # parsed_args
                nothing, # conservation_checks
                (Int(0), Int(1000)), # tspan
                Int(100), # t
                Int(100), # Δt_cpl
                (;), # surface_masks
                (;), # model_sims
                (;), # mode
                (dg_2d,),
                (;), # callbacks
                (;), # dirs
                nothing, # turbulent_fluxes
                nothing, # thermo_params
            )
            Diagnostics.accumulate_diagnostics!(cs)
            @test cs.diagnostics[1].field_vector[1] == expected_results[c_i][1]
            Diagnostics.accumulate_diagnostics!(cs)
            Diagnostics.pre_save(cs.diagnostics[1].operations.accumulate, cs, cs.diagnostics[1])
            @test cs.diagnostics[1].field_vector[1] == expected_results[c_i][2]

            Diagnostics.post_save(cs.diagnostics[1].operations.accumulate, cs, cs.diagnostics[1])
            @test cs.diagnostics[1].field_vector[1] == expected_results[c_i][3]
        end
    end
end
