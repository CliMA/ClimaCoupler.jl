export CoupledSimulation, step!, run!

"""
    CoupledSimulation

An abstract type representing a coupled simulation.
"""
abstract type CoupledSimulation end

"""
    run!(::CoupledSimulation)

A simple outer timestepping loop for coupled system runs.

This will be formalized when the run! functionality for component
models is implemented so to have a consistent interface.
"""
function run!(sim::CoupledSimulation)
    clock = sim.clock
    while !stop_time_exceeded(clock)
        step!(sim, clock.dt)
        tick!(clock)
    end
end

"""
    step!(sim, dt)

Advances a simulation by `dt`.

Note that `dt` is not necessarily the simulation's timestep length;
a simuation could take several shorter steps that total to `dt`.
"""
function step!(sim, dt) end




#=
__________________

function run!(::Sequential, coupler_sim)

    for t in timesteps

        cpl_bdry_flux = coupler_callback(coupler_sim) # ref

        for component in coupler_sim
            pull!(component.boundary, cpl_bdry_flux)
        end
        for
            step!
        end
        for
            push!(cpl_state, component.state) # if needed, but not needed if no accumulation (since no storage in coupler)
        end
    end
end

function run!(::Sequential, coupler_sim)

    driver_component = coupler_sim.component[1]
    for t in timesteps

        calc_fluxes!(driver_component.p, coupler_sim) # calc fluxes and save in atmos.p at every coupler_sim dt

        for component in coupler_sim
            pull!(component.boundary, driver_component.p)
            step!
            push!(driver_component.state, component.state) # if needed, but not needed if no accumulation (since no storage in coupler)
        end
    end
end

function calc_fluxes!(driver_component.p, coupler_sim)
    driver_component.p.flux = calc_fluxes(coupler_sim)
end

______

calc_fluxes!(driver_component.p, coupler_sim) # calc fluxes and save in atmos.p

driver_sim = driver_sim_init(..., callback = calc_fluxes!) # flux calculation at every driver_sim dt

coupler_sim = (driver_sim, sim2, sim3, ...)

function run!(::Sequential, coupler_sim)

    driver_component = coupler_sim.component[1]
    for t in timesteps

        for component in coupler_sim
            pull!(component.boundary, driver_component.p)
            step!
            push!(driver_component.state, component.state) # if needed, but not needed if no accumulation (since no storage in coupler)
        end
    end
end

function calc_fluxes!(driver_component.p, coupler_sim)
    driver_component.p.flux = calc_fluxes(coupler_sim)
end



Set up coupler:
1. Registration phase - tell the coupler what fields
     are available to other models via `coupler_get`
...
2. Driver model gets some fields "from" the coupler, calculates flux and
    stores in some flux variable it owns
3. Other models using that flux can access it via `coupler_get` which
    reaches into driver model and regrids if needed

=#