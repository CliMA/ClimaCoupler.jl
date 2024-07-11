#=
    Unit tests for ClimaCoupler Utilities module
=#

import Dates
using Test
import ClimaComms
using ClimaCoupler: Utilities, TestHelper
using ClimaCore: Fields

for FT in (Float32, Float64)
    @testset "test swap_space!" begin
        space1 = TestHelper.create_space(FT, R = FT(6371e3))
        space2 = TestHelper.create_space(FT, R = FT(6371e3))

        field1 = ones(space1)
        field2 = ones(space2)

        field2 = Utilities.swap_space!(space2, field1)

        @test parent(field1) == parent(field2)
        @test axes(field2) == space2
    end

    @testset "test current_date" begin
        date0 = date = Dates.DateTime("19790321", Dates.dateformat"yyyymmdd")
        dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])
        tspan = (Int(1), Int(90 * 86400)) # Jan-Mar
        Δt_cpl = 1 * 24 * 3600

        # Fill in only the necessary parts of the simulation
        cs = Utilities.CoupledSimulation{FT}(
            ClimaComms.SingletonCommsContext(), # comms_ctx
            dates, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            tspan, # tspan
            Int(0), # t
            Int(Δt_cpl), # Δt_cpl
            (;), # surface_masks
            (;), # model_sims
            (;), # mode
            (), # diagnostics
        )

        for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])
            @test Utilities.current_date(cs, t) == date0 + Dates.Second(t)
        end
    end
end
