import ClimaAnalysis
import Dates
import ClimaCalibrate
import GeoMakie
import Makie

import ClimaCalibrate: EnsembleBuilder

include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/observation_utils.jl"))

"""
    ClimaCalibrate.observation_map(iteration)

Return G ensemble for an `iteration`.

G ensemble represents the concatenated forward model evaluations from all
ensemble members, arranged horizontally. Each individual forward model
evaluation corresponds to preprocessed, flattened simulation data from a single
ensemble member that has been matched to the corresponding observational data.
"""
function ClimaCalibrate.observation_map(iteration)
    output_dir = CALIBRATE_CONFIG.output_dir
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))

    g_ens_builder = EnsembleBuilder.GEnsembleBuilder(ekp)

    G_ensemble = Array{Float64}(undef, single_obs_len, ensemble_size)
    for m = 1:ensemble_size
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "global_diagnostics/output_active")
        @info "Processing member $m: $simdir_path"
        try
                process_member_data(g_ens_builder, simdir_path, col_idx)
        catch e
            @error "Error processing member $m, filling observation map entry with NaNs" exception =
                e
            G_ensemble[:, m] .= NaN
        end
    end

    return G_ensemble
end

"""
    process_member_data(diagnostics_folder_path, short_names, current_minibatch)

Process the data of a single ensemble member and return a single column of the
G ensemble matrix.
"""
function process_member_data(g_ens_builder, diagnostics_folder_path, col_idx)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    nelements = CALIBRATE_CONFIG.nelements
    @info "Short names: $short_names"

    # For now, hard code the data to daily data since it seems the simplest
    # TODO: Fetch daily data from simdir
    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in short_names
        # This assumes that there are not multiple reductions for the same
        # variable!
        var = get(simdir, short_name)
        var = preprocess_var(var, first(first(sample_date_ranges)))
        EnsembleBuilder.fill_g_ens_col!(g_ens_builder, col_idx, var)
    end

    return nothing
end

"""
    preprocess_var(var::ClimaAnalysis.OutputVar, reference_date)

Preprocess `var` before flattening for G ensemble matrix.

For "pr", weekly sums are computed. For "tas" and "mslp", weekly means are
computed from daily means. The daily means are computing starting from
`reference_date`.

This function assumes that the data is daily.
"""
function preprocess_var(var::ClimaAnalysis.OutputVar, reference_date)
    # TODO: Check for sign of pr
    # TODO: Check for units of everything
    var = shift_to_previous_week(var)
    if ClimaAnalysis.short_name(var) == "pr"
        # TODO: Check that the sign are the same as the observational data
        return compute_weekly_sum_from_daily_mean(var, reference_date)
    elseif ClimaAnalysis.short_name(var) == "tas"
        return compute_weekly_mean_from_daily_mean(var, reference_date)
    elseif ClimaAnalysis.short_name(var) == "mslp"
        return compute_weekly_mean_from_daily_mean(var, reference_date)
    else
        error("Do not know how to preprocess OutputVar with short name $(ClimaAnalysis.short_name(var))")
    end
end

"""
    ClimaCalibrate.analyze_iteration(ekp,
                                     g_ensemble,
                                     prior,
                                     output_dir,
                                     iteration)

Analyze an iteration by plotting the bias plots, constrained parameters over
iterations, and errors over iterations and time.
"""
function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
    plot_constrained_params_and_errors(plot_output_path, ekp, prior)

    # Plot ERA5 bias plots for only the first ensemble member
    # This can take a while to plot, so we plot only one of the members.
    # We choose the first ensemble member because the parameters for the first
    # ensemble member are supposed to be the mean of the parameters of the
    # ensemble members if it is EKP.TransformUnscented
    # TODO: Add plotting here for the first ensemble member!
end

"""
    plot_constrained_params_and_errors(output_dir, ekp, prior)

Plot the constrained parameters and errors from `ekp` and `prior` and save
them to `output_dir`.
"""
function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i = 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size+1], ekp, error_metric = "loss")
    EKP.Visualize.plot_error_over_time(fig[1, dim_size+2], ekp, error_metric = "loss")
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
end
