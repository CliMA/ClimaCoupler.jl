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

    for m = 1:EKP.get_N_ens(ekp)
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "wxquest_diagedmf/output_active")
        @info "Processing member $m: $simdir_path"
        process_member_data(g_ens_builder, simdir_path, m, iteration)
    end

    return 
end

"""
    process_member_data(diagnostics_folder_path, short_names, current_minibatch)

Process the data of a single ensemble member and return a single column of the
G ensemble matrix.
"""
function process_member_data(g_ens_builder, diagnostics_folder_path, col_idx, iteration)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    @info "Short names: $short_names"

    # For now, hard code the data to daily data since it seems the simplest
    # TODO: Fetch daily data from simdir
    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in short_names
        # This assumes that there are not multiple reductions for the same
        # variable!
        short_name = short_name * "_1week"
        short_name == "tas_1week" && (short_name = "ta_1week")
        var = get(simdir, short_name)
        var.attributes["short_name"] = replace(var.attributes["short_name"], "_1week" => "")
        # var = preprocess_var(var, first(sample_date_ranges[iteration+1]))
        var = preprocess_var(var, DateTime(2024,9,29))
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
    @show var |> dates
    # @assert ClimaAnalysis.short_name(var) in CALIBRATE_CONFIG.short_names
    if ClimaAnalysis.short_name(var) == "pr"
        # Turn weekly average into sum
        var.data .*= 7
    end
    @show reference_date
    var = window(var, "time"; left = reference_date, right = reference_date)
    @show var |> dates
    return var
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