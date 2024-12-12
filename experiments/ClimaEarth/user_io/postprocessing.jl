## helpers for user-specified IO
include("debug_plots.jl")
include("diagnostics_plots.jl")

"""
    postprocess_sim(sim_mode::AbstractSlabplanetSimulationMode, cs, postprocessing_vars)

If they exist, perform conservation checks, for any slabplanet simulation, and then call
`common_postprocessing` to perform common postprocessing tasks that are common to all simulation types.
"""
function postprocess_sim(sim_mode::AbstractSlabplanetSimulationMode, cs, postprocessing_vars)
    (; conservation_softfail,) = postprocessing_vars

    common_postprocessing(cs, postprocessing_vars)

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
end

"""
    postprocess_sim(sim_mode::AMIPMode, cs, postprocessing_vars)

Conditionally plot AMIP diagnostics and call `common_postprocessing` to perform
postprocessing tasks that are common to all simulation types.
"""
function postprocess_sim(sim_mode::AMIPMode, cs, postprocessing_vars)
    (; use_coupler_diagnostics, output_default_diagnostics, t_end) = postprocessing_vars

    common_postprocessing(cs, postprocessing_vars)

    if use_coupler_diagnostics
        ## plot data that correspond to the model's last save_hdf5 call (i.e., last month)
        @info "AMIP plots"

        ## ClimaESM
        include("user_io/diagnostics_plots.jl")

        # define variable names and output directories for each diagnostic
        amip_short_names_atmos = ["ta", "ua", "hus", "clw", "pr", "ts", "toa_fluxes_net"]
        amip_short_names_coupler = ["F_turb_energy"]
        output_dir_coupler = cs.dirs.output

        # Check if all output variables are available in the specified directories
        make_diagnostics_plots(
            atmos_output_dir,
            cs.dirs.artifacts,
            short_names = amip_short_names_atmos,
            output_prefix = "atmos_",
        )
        make_diagnostics_plots(
            output_dir_coupler,
            cs.dirs.artifacts,
            short_names = amip_short_names_coupler,
            output_prefix = "coupler_",
        )
    end

    # Check this because we only want monthly data for making plots
    if t_end > 84600 * 31 * 3 && output_default_diagnostics
        include("leaderboard/leaderboard.jl")
        leaderboard_base_path = cs.dirs.artifacts
        compute_leaderboard(leaderboard_base_path, atmos_output_dir)
        compute_pfull_leaderboard(leaderboard_base_path, atmos_output_dir)
    end

    !isnothing(cs.amip_diags_handler) &&
        map(diag -> close(diag.output_writer), cs.amip_diags_handler.scheduled_diagnostics)
end

"""
    common_postprocessing(cs, postprocessing_vars)

Perform postprocessing common to all simulation types.
"""
function common_postprocessing(cs, postprocessing_vars)
    (; plot_diagnostics, atmos_output_dir) = postprocessing_vars
    if plot_diagnostics
        @info "Plotting diagnostics"
        include("user_io/diagnostics_plots.jl")
        make_diagnostics_plots(atmos_output_dir, cs.dirs.artifacts)
    end

    # plot all model states and coupler fields (useful for debugging)
    !CA.is_distributed(cs.comms_ctx) && debug(cs, cs.dirs.artifacts)
end
