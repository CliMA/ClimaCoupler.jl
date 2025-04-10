## helpers for user-specified IO
include("debug_plots.jl")
include("diagnostics_plots.jl")
include("../leaderboard/leaderboard.jl")

"""
    postprocess_sim(cs, postprocessing_vars)

Perform all postprocessing operations. This includes plotting all available
diagnostics, plotting all model states and coupler fields for debugging,
producing the leaderboard if monthly data is available, performing
conservation checks if enabled, and closing all diagnostics file writers.
"""
function postprocess_sim(cs, postprocessing_vars)
    (; conservation_softfail) = postprocessing_vars
    output_dir = cs.dirs.output
    artifact_dir = cs.dirs.artifacts
    coupler_output_dir = joinpath(output_dir, "coupler")
    atmos_output_dir = joinpath(output_dir, "clima_atmos")
    land_output_dir = joinpath(output_dir, "clima_land")

    # Plot generic diagnostics
    @info "Plotting diagnostics for coupler, atmos, and land"
    make_diagnostics_plots(coupler_output_dir, artifact_dir, output_prefix = "coupler_")
    make_diagnostics_plots(atmos_output_dir, artifact_dir, output_prefix = "atmos_")
    make_diagnostics_plots(land_output_dir, artifact_dir, output_prefix = "land_")

    # Plot all model states and coupler fields (useful for debugging)
    # TODO can't make debug plots because we have some NaNs
    # !CA.is_distributed(cs.comms_ctx) && debug(cs, artifact_dir)

    # If we have enough data (in time, but also enough variables), plot the leaderboard.
    # We need pressure to compute the leaderboard.
    pressure_in_output = "pfull" in CAN.available_vars(CAN.SimDir(atmos_output_dir))
    if pressure_in_output
        times = CAN.times(get(CAN.SimDir(atmos_output_dir), "pfull"))
        t_end = times[end]
        if t_end > 84600 * 31 * 3 # 3 months for spin up
            leaderboard_base_path = artifact_dir
            compute_leaderboard(leaderboard_base_path, atmos_output_dir)
            compute_pfull_leaderboard(leaderboard_base_path, atmos_output_dir)
        end
    end

    # Perform conservation checks if they exist
    if !isnothing(cs.conservation_checks)
        @info "Conservation Check Plots"
        plot_global_conservation(
            cs.conservation_checks.energy,
            cs,
            conservation_softfail,
            figname1 = joinpath(cs.dirs.artifacts, "total_energy_bucket.png"),
            figname2 = joinpath(cs.dirs.artifacts, "total_energy_log_bucket.png"),
        )
        plot_global_conservation(
            cs.conservation_checks.water,
            cs,
            conservation_softfail,
            figname1 = joinpath(cs.dirs.artifacts, "total_water_bucket.png"),
            figname2 = joinpath(cs.dirs.artifacts, "total_water_log_bucket.png"),
        )
    end

    # Close all diagnostics file writers
    !isnothing(cs.diags_handler) && map(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)
end
