import ClimaAtmos
import ClimaComms
import ClimaCoupler: Checkpointer, CoupledSimulation, Input, Interfacer, step!
using Test

ClimaComms.@import_required_backends

slab_value(field) = only(Array(parent(field)))

function slab_heat_fluxes(slab)
    (; SW_d, LW_d, F_turb_energy) = slab.integrator.p
    return map(slab_value, (; SW_d, LW_d, F_turb_energy))
end

function slab_window_heat(slab, forcing, start_temperature, end_temperature)
    (; h, ρ, c, α, ϵ, σ) = slab.params
    (; SW_d, LW_d, F_turb_energy) = forcing
    dt = Float64(slab.integrator.dt)
    heat_capacity = h * ρ * c
    input = (1 - α) * SW_d + ϵ * LW_d - F_turb_energy
    emission_start = ϵ * σ * start_temperature^4
    emission_end = ϵ * σ * end_temperature^4
    # Estimate the midpoint temperature from the two endpoint rates, then integrate emission.
    midpoint =
        (start_temperature + end_temperature) / 2 +
        dt * (emission_end - emission_start) / (8 * heat_capacity)
    mean_emission = (emission_start + 4 * ϵ * σ * midpoint^4 + emission_end) / 6
    return dt * (input - mean_emission)
end

function finish_slab_run!(cs)
    slab = cs.model_sims.ocean_sim
    temperatures = Float64[]
    window_heat = Float64[]
    (; h, ρ, c) = slab.params
    initial_temperature = slab_value(slab.integrator.u.T_sfc)
    try
        while cs.t[] < cs.tspan[end]
            clock_before = slab.integrator.t
            temperature_before = slab_value(slab.integrator.u.T_sfc)
            # Read the forcing before the coupling step refreshes the cache.
            forcing = slab_heat_fluxes(slab)
            step!(cs)
            temperature = slab_value(slab.integrator.u.T_sfc)
            if slab.integrator.t != clock_before
                push!(
                    window_heat,
                    slab_window_heat(slab, forcing, temperature_before, temperature),
                )
            end
            push!(temperatures, temperature)
        end
    finally
        foreach(Interfacer.close_output_writers, cs.model_sims)
    end
    heat_storage = h * ρ * c * (last(temperatures) - initial_temperature)
    return (; temperatures, window_heat, heat_storage)
end

@testset "Restart inside a slow slab-ocean flux window at $checkpoint_s s" for (
    checkpoint_s,
    flux_samples,
) in (
    (480, 2),
    (600, 0),
)
    mktempdir() do dir
        config_file = joinpath(@__DIR__, "../../config/ci_configs/scm_slabplanet_aqua.yml")
        config = Input.get_coupler_config_dict(config_file)
        config["FLOAT_TYPE"] = "Float64"
        config["ocean_model"] = "slab"
        config["coupler_output_dir"] = dir
        config["job_id"] = "slow_slab_reference"
        config["enable_diagnostics"] = false
        config["use_coupler_diagnostics"] = false
        config["print_config_dict"] = false
        config["t_end"] = "1080secs"
        config["dt_atmos"] = config["dt_land"] = config["dt_seaice"] = config["dt_cpl"]
        config["dt_ocean"] = "360secs"
        delete!(config, "dt")
        Input.parse_component_dts!(config)

        reference = CoupledSimulation(config)
        for _ in 1:div(checkpoint_s, 120)
            step!(reference)
        end
        slab = reference.model_sims.ocean_sim
        @test reference.flux_accumulators.ocean_sim.n_steps[] == flux_samples
        checkpoint_fluxes = slab_heat_fluxes(slab)
        Checkpointer.checkpoint_sims(reference)
        expected = finish_slab_run!(reference)

        config["job_id"] = "slow_slab_restart"
        config["restart_dir"] = reference.dir_paths.checkpoints_dir
        config["restart_t"] = checkpoint_s
        restarted = CoupledSimulation(config)
        restarted_slab = restarted.model_sims.ocean_sim
        if flux_samples > 0
            @test Float64(restarted_slab.integrator.t) == 360
        else
            @test slab_heat_fluxes(restarted_slab) == checkpoint_fluxes
        end
        actual = finish_slab_run!(restarted)
        @test actual.temperatures == expected.temperatures
        @test actual.window_heat == expected.window_heat
        @test abs(sum(actual.window_heat) - actual.heat_storage) < 1e-4
    end
end
