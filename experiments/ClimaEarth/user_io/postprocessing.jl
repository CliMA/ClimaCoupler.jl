## helpers for user-specified IO
include("debug_plots.jl")
include("diagnostics_plots.jl")
include("../leaderboard/leaderboard.jl")
include("../leaderboard/test_rmses.jl")

"""
    postprocess_sim(cs, postprocessing_vars)

Perform all postprocessing operations. This includes plotting all available
diagnostics, plotting all model states and coupler fields for debugging,
producing the leaderboard if monthly data is available, performing
conservation checks if enabled, and closing all diagnostics file writers.
"""
function postprocess_sim(cs, postprocessing_vars)
    (; conservation_softfail, rmse_check) = postprocessing_vars
    (;
        coupler_output_dir,
        atmos_output_dir,
        land_output_dir,
        ocean_output_dir,
        artifacts_dir,
    ) = cs.dir_paths

    # Plot generic diagnostics
    @info "Plotting diagnostics for coupler, atmos, land, and ocean"
    make_diagnostics_plots(coupler_output_dir, artifacts_dir, output_prefix = "coupler_")
    make_diagnostics_plots(atmos_output_dir, artifacts_dir, output_prefix = "atmos_")
    make_diagnostics_plots(land_output_dir, artifacts_dir, output_prefix = "land_")
    make_ocean_diagnostics_plots(ocean_output_dir, artifacts_dir, output_prefix = "ocean_")

    # Plot all model states and coupler fields (useful for debugging)
    ClimaComms.context(cs) isa ClimaComms.SingletonCommsContext && debug(cs, artifacts_dir)

    # If we have enough data (in time, but also enough variables), plot the leaderboard.
    # We need pressure to compute the leaderboard.
    simdir = CAN.SimDir(atmos_output_dir)
    if !isempty(simdir)
        pressure_in_output = "pfull" in CAN.available_vars(simdir)
        times = CAN.times(get(simdir, first(CAN.available_vars(simdir))))
        t_end = times[end]
        if t_end > 84600 * 31 * 3 # 3 months for spin up
            leaderboard_base_path = artifacts_dir
            compute_leaderboard(leaderboard_base_path, atmos_output_dir, 3)
            rmse_check && test_rmse_thresholds(atmos_output_dir, 3)
            pressure_in_output &&
                compute_pfull_leaderboard(leaderboard_base_path, atmos_output_dir, 3)
        end
    end

    # Perform conservation checks if they exist
    if !isnothing(cs.conservation_checks)
        @info "Conservation Check Plots"
        plot_global_conservation(
            cs.conservation_checks.energy,
            cs,
            conservation_softfail,
            figname1 = joinpath(artifacts_dir, "total_energy_bucket.png"),
            figname2 = joinpath(artifacts_dir, "total_energy_log_bucket.png"),
        )
        plot_global_conservation(
            cs.conservation_checks.water,
            cs,
            conservation_softfail,
            figname1 = joinpath(artifacts_dir, "total_water_bucket.png"),
            figname2 = joinpath(artifacts_dir, "total_water_log_bucket.png"),
        )
    end

    # Close all diagnostics file writers
    !isnothing(cs.diags_handler) &&
        map(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)
end

"""
    simulated_years_per_day(cs, walltime)

Compute the simulated years per walltime day for the given coupled simulation `cs`, assuming
that the simulation took `walltime`.
"""
function simulated_years_per_day(cs, walltime)
    simulated_seconds_per_second = float(cs.tspan[end] - cs.tspan[begin]) / walltime
    return simulated_seconds_per_second / 365.25
end

"""
    walltime_per_coupling_step(cs, walltime)

Compute the average walltime needed to take one step for the given coupled simulation `cs`,
assuming that the simulation took `walltime`. The result is in seconds.
"""
function walltime_per_coupling_step(cs, walltime)
    n_coupling_steps = (cs.tspan[end] - cs.tspan[begin]) / cs.Δt_cpl
    return walltime / n_coupling_steps
end


"""
    save_sypd_walltime_to_disk(cs, walltime)

Save the computed `sypd`, `walltime_per_coupling_step`,
and memory usage to text files in the `artifacts` directory.
"""
function save_sypd_walltime_to_disk(cs, walltime)
    if ClimaComms.iamroot(ClimaComms.context(cs))
        sypd = simulated_years_per_day(cs, walltime)
        walltime_per_step = walltime_per_coupling_step(cs, walltime)

        open(joinpath(cs.dir_paths.artifacts_dir, "sypd.txt"), "w") do sypd_filename
            println(sypd_filename, "$sypd")
        end

        open(
            joinpath(cs.dir_paths.artifacts_dir, "walltime_per_step.txt"),
            "w",
        ) do walltime_per_step_filename
            println(walltime_per_step_filename, "$(walltime_per_step)")
        end
    end
    return nothing
end
