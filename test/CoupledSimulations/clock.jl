using Test
using CouplerMachine: Clock, tick!, stop_time_exceeded

@testset "Clock" begin
    time_info = (start = 0.0, dt = 0.5, stop = 2.0)
    clock = Clock(time_info...)

    tick!(clock)
    @test clock.time == time_info.dt
    while !stop_time_exceeded(clock)
        tick!(clock)
    end
    @test clock.time == time_info.stop
end
