#=
    Unit tests for ClimaCoupler TimeManager module
=#

using Test
using Dates
using ClimaCoupler: Interfacer, TimeManager
using ClimaComms

for FT in (Float32, Float64)
    @testset "test current_date" begin
        date0 = date = DateTime("19790321", dateformat"yyyymmdd")
        dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])
        tspan = (Int(1), Int(90 * 86400)) # Jan-Mar
        Δt_cpl = 1 * 24 * 3600

        # Fill in only the necessary parts of the simulation
        cs = Interfacer.CoupledSimulation{FT}(
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
            @test TimeManager.current_date(cs, t) == date0 + Dates.Second(t)
        end
    end
end

@testset "test strdate_to_datetime" begin
    @test TimeManager.strdate_to_datetime("19000101") == Dates.DateTime(1900, 1, 1)
    @test TimeManager.strdate_to_datetime("00000101") == Dates.DateTime(0, 1, 1)
end

@testset "test datetime_to_strdate" begin
    @test TimeManager.datetime_to_strdate(Dates.DateTime(1900, 1, 1)) == "19000101"
    @test TimeManager.datetime_to_strdate(Dates.DateTime(0, 1, 1)) == "00000101"
end

@testset "trigger_callback" begin
    date0 = date = DateTime("19790321", dateformat"yyyymmdd")
    dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])

    cs = Interfacer.CoupledSimulation{Float64}(
        nothing, # comms_ctx
        dates, # dates
        nothing, # boundary_space
        nothing, # fields
        nothing, # parsed_args
        nothing, # conservation_checks
        (Int(0), Int(1000)), # tspan
        Int(200), # t
        Int(200), # Δt_cpl
        (;), # surface_masks
        (;), # model_sims
        (;), # mode
        (), # diagnostics
    )
    @test TimeManager.trigger_callback(cs, TimeManager.Monthly()) == true
end
