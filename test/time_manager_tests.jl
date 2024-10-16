#=
    Unit tests for ClimaCoupler TimeManager module
=#
import Test: @testset, @test
import Dates
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCoupler: Interfacer, TimeManager

for FT in (Float32, Float64)
    @testset "test current_date" begin
        date0 = date = Dates.DateTime("19790321", Dates.dateformat"yyyymmdd")
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
            (;), # model_sims
            (;), # mode
            (;), # callbacks
            (;), # dirs
            nothing, # turbulent_fluxes
            nothing, # thermo_params
            nothing, # amip_diags_handler
        )

        for t in ((tspan[1]+Δt_cpl):Δt_cpl:tspan[end])
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
    FT = Float64
    date0 = date = Dates.DateTime("19790321", Dates.dateformat"yyyymmdd")
    dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])

    cs = Interfacer.CoupledSimulation{FT}(
        nothing, # comms_ctx
        dates, # dates
        nothing, # boundary_space
        nothing, # fields
        nothing, # parsed_args
        nothing, # conservation_checks
        (Int(0), Int(1000)), # tspan
        Int(200), # t
        Int(200), # Δt_cpl
        (;), # model_sims
        (;), # mode
        (;), # callbacks
        (;), # dirs
        nothing, # turbulent_fluxes
        nothing, # thermo_params
        nothing, # amip_diags_handler
    )
    @test TimeManager.trigger_callback(cs, TimeManager.Monthly()) == true
end

@testset "trigger_callback!" begin
    FT = Float64
    date0 = date = Dates.DateTime("19790321", Dates.dateformat"yyyymmdd")
    dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])

    function counter_func(cs, cb)
        cb.data .+= 1
    end
    twhohourly_inactive = TimeManager.HourlyCallback{FT}(dt = 2, ref_date = [date0])
    twhohourly_nothing = TimeManager.HourlyCallback{FT}(dt = 2, ref_date = [date0], active = true)
    twhohourly_counter =
        TimeManager.HourlyCallback{FT}(dt = 2, ref_date = [date0], func = counter_func, data = [0], active = true)
    monthly_counter =
        TimeManager.MonthlyCallback{FT}(func = counter_func, ref_date = [date0], data = [0], active = true)

    cs = Interfacer.CoupledSimulation{FT}(
        nothing, # comms_ctx
        dates, # dates
        nothing, # boundary_space
        nothing, # fields
        nothing, # parsed_args
        nothing, # conservation_checks
        (Int(0), Int(1000)), # tspan
        Int(200), # t
        Int(200), # Δt_cpl
        (;), # model_sims
        (;), # mode
        (;
            twhohourly_inactive = twhohourly_inactive,
            twhohourly_nothing = twhohourly_nothing,
            twhohourly_counter = twhohourly_counter,
            monthly_counter = monthly_counter,
        ), # callbacks
        (;), # dirs
        nothing, # turbulent_fluxes
        nothing, # thermo_params
        nothing, # amip_diags_handler
    )

    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_inactive)
    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_nothing)
    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_counter)
    TimeManager.trigger_callback!(cs, cs.callbacks.monthly_counter)
    @test cs.callbacks.twhohourly_inactive.ref_date[1] == date0
    @test cs.callbacks.twhohourly_nothing.ref_date[1] == cs.dates.date0[1] + Dates.Hour(2)
    @test cs.callbacks.twhohourly_counter.ref_date[1] == cs.dates.date0[1] + Dates.Hour(2)
    @test cs.callbacks.monthly_counter.ref_date[1] == cs.dates.date0[1] + Dates.Month(1)
    @test cs.callbacks.twhohourly_counter.data[1] == 1
    @test cs.callbacks.monthly_counter.data[1] == 1

    cs.dates.date .+= Dates.Hour(2)
    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_inactive)
    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_nothing)
    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_counter)
    TimeManager.trigger_callback!(cs, cs.callbacks.monthly_counter)
    @test cs.callbacks.twhohourly_inactive.ref_date[1] == date0
    @test cs.callbacks.twhohourly_nothing.ref_date[1] == cs.dates.date0[1] + Dates.Hour(4)
    @test cs.callbacks.twhohourly_counter.ref_date[1] == cs.dates.date0[1] + Dates.Hour(4)
    @test cs.callbacks.monthly_counter.ref_date[1] == cs.dates.date0[1] + Dates.Month(1)
    @test cs.callbacks.twhohourly_counter.data[1] == 2
    @test cs.callbacks.monthly_counter.data[1] == 1

    cs.dates.date .+= Dates.Month(1)
    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_inactive)
    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_nothing)
    TimeManager.trigger_callback!(cs, cs.callbacks.twhohourly_counter)
    TimeManager.trigger_callback!(cs, cs.callbacks.monthly_counter)
    @test cs.callbacks.twhohourly_inactive.ref_date[1] == date0
    @test cs.callbacks.twhohourly_inactive.ref_date[1] == date0
    @test cs.callbacks.twhohourly_nothing.ref_date[1] == cs.dates.date0[1] + Dates.Hour(6)
    @test cs.callbacks.twhohourly_counter.ref_date[1] == cs.dates.date0[1] + Dates.Hour(6)
    @test cs.callbacks.monthly_counter.ref_date[1] == cs.dates.date0[1] + Dates.Month(2)
    @test cs.callbacks.twhohourly_counter.data[1] == 3
    @test cs.callbacks.monthly_counter.data[1] == 2

end

# TimeManager
@testset "update_firstdayofmonth!" begin
    FT = Float64
    date0 = date = Dates.DateTime("19790321", Dates.dateformat"yyyymmdd")
    dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])

    cs = Interfacer.CoupledSimulation{FT}(
        nothing, # comms_ctx
        dates, # dates
        nothing, # boundary_space
        nothing, # fields
        nothing, # parsed_args
        nothing, # conservation_checks
        (Int(0), Int(1000)), # tspan
        Int(200), # t
        Int(200), # Δt_cpl
        (;), # model_sims
        (;), # mode
        (;), # callbacks
        (;), # dirs
        nothing, # turbulent_fluxes
        nothing, # thermo_params
        nothing, # amip_diags_handler
    )

    TimeManager.update_firstdayofmonth!(cs, nothing)
    @test cs.dates.date1[1] == Dates.firstdayofmonth(date0) + Dates.Month(1)

end
