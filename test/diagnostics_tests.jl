#= 
    Unit tests for ClimaCoupler Diagnostics module
=#
using Test
using Dates
using ClimaComms
using ClimaCoupler: Utilities, TestHelper

FT = Float64
get_var(cs::CoupledSimulation, ::Val{:x}) = FT(1)

@testset "init_diagnostics" begin
    names = (:x, :y)
    space = create_space(FT)
    dg = init_diagnostics(names, space)
    @test typeof(dg) == DiagnosticsGroup{EveryTimestep}
end

@testset "accumulate_diagnostics!, collect_diags, iterate_operations, operation{accumulation{TimeMean, Nothing}}" begin
    cases = (nothing, TimeMean([Int(0)]))
    expected_results = (FT(2), FT(3))
    for (c_i, case) in enumerate(cases)
        names = (:x,)
        space = create_space(FT)
        dg_2d = init_diagnostics(names, space, save = EveryTimestep(), operations = (; accumulate = case))
        dg_2d.field_vector .= FT(2)
        cs = Utilities.CoupledSimulation(
            nothing, # comms_ctx
            nothing, # tspan
            nothing, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            FT(100), # t
            FT(100), # Δt_cpl
            (;), # surface_masks
            (;), # model_sims
            (;), # mode
            (dg_2d,),
        )
        accumulate_diagnostics!(cs)
        @test cs.diagnostics[1].field_vector[1] == expected_results[c_i]

    end
end

@testset "save_diagnostics" begin
    test_dir = "diag_test_dir"
    nmes = (:x,)
    space = create_space(FT)
    dg_2d = init_diagnostics(
        nmes,
        space,
        save = EveryTimestep(),
        operations = (; accumulate = TimeMean([Int(0)])),
        output_dir = test_dir,
    ) # or use accumulate = nothing for snapshop save
    cs = Utilities.CoupledSimulation(
        ClimaComms.SingletonCommsContext(), # comms_ctx
        nothing,# tspan
        (date = [DateTime(0, 2)], date1 = [DateTime(0, 1)]), # dates
        nothing, # boundary_space
        nothing, # fields
        nothing, # parsed_args
        nothing, # conservation_checks
        FT(100), # t
        FT(100), # Δt_cpl
        (;), # surface_masks
        (;), # model_sims
        (;), # mode
        (dg_2d,), # diagnostics
    )
    save_diagnostics(cs, cs.diagnostics[1])
    file = filter(x -> endswith(x, ".hdf5"), readdir(test_dir))
    @test !isempty(file)
    rm(test_dir; recursive = true, force = true)

end

@testset "save_time_format" begin
    date = DateTime(1970, 2, 1, 0, 1)
    unix = save_time_format(date, Monthly())
    @test unix == 0
end

@testset "pre_save{TimeMean, Nothing}, post_save" begin
    cases = (nothing, TimeMean([Int(0)]))
    expected_results = ((FT(3), FT(1), FT(1)), (FT(4), FT(2.5), FT(0)))

    for (c_i, case) in enumerate(cases)
        names = (:x,)
        space = create_space(FT)
        dg_2d = init_diagnostics(nmes, space, save = EveryTimestep(), operations = (; accumulate = case))
        dg_2d.field_vector .= FT(3)
        cs = Utilities.CoupledSimulation(
            nothing, # comms_ctx
            nothing, # tspan
            nothing, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            FT(100), # t
            FT(100), # Δt_cpl
            (;), # surface_masks
            (;), # model_sims
            (;), # mode
            (dg_2d,),
        )
        accumulate_diagnostics!(cs)
        @test cs.diagnostics[1].field_vector[1] == expected_results[c_i][1]
        accumulate_diagnostics!(cs)
        pre_save(cs.diagnostics[1].operations.accumulate, cs, cs.diagnostics[1])
        @test cs.diagnostics[1].field_vector[1] == expected_results[c_i][2]

        post_save(cs.diagnostics[1].operations.accumulate, cs, cs.diagnostics[1])
        @test cs.diagnostics[1].field_vector[1] == expected_results[c_i][3]
    end
end
