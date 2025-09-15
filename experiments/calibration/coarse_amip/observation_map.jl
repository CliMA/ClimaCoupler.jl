import ClimaAnalysis
import Dates
import ClimaCalibrate
import GeoMakie
import Makie
import ClimaAnalysis.Utils: kwargs as ca_kwargs

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
        try
            process_member_data!(g_ens_builder, simdir_path, m, iteration)
        catch e
            @error "Ensemble member $m failed" exception =
            (e, catch_backtrace())
        end
    end

    return g_ens_builder.g_ens
end

"""
    process_member_data!(diagnostics_folder_path, short_names, current_minibatch)

Process the data of a single ensemble member and return a single column of the
G ensemble matrix.
"""
function process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx, iteration)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    @info "Short names: $short_names"

    # For now, hard code the data to daily data since it seems the simplest
    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in short_names
        # This assumes that there are not multiple reductions for the same
        # variable!
        if short_name == "tas"
            short_name = "ts"
        end
        short_name = short_name * "_1week" # TODO: Don't hardcode this
        var = get(simdir, short_name)
        var.attributes["short_name"] = replace(short_name, "_1week" => "")

        var = preprocess_var(var, first(sample_date_ranges[iteration+1]))

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
    # @assert ClimaAnalysis.short_name(var) in CALIBRATE_CONFIG.short_names
    if ClimaAnalysis.short_name(var) == "pr"
        # Turn weekly average into sum
        var.data .*= 7
    elseif ClimaAnalysis.short_name(var) == "ts"
        var.attributes["short_name"] = "tas"
    end
    
    var = window(var, "time"; left = reference_date, right = reference_date)
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

function plot_variables(simdir)
    pr = get(simdir, "pr_1week")
    mslp = get(simdir, "mslp_1week")
    ts = get(simdir, "ts_1week")

    pr.data .*= 7
    pr.attributes["long_name"] = "Precipitation, accumulation within 1 Week"

    fig = GeoMakie.Figure(size = (1000, 3*500))

    for (i, var) in enumerate((pr, mslp, ts))
        var = slice(var, time = last(dates(var)))
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[i, 1], var; more_kwargs = Dict(:plot => ca_kwargs(colormap = :viridis)))
    end

    GeoMakie.save(joinpath(simdir.simulation_path, "vars.png"), fig)
end

function plot_bias(simdir)
    pr = get(simdir, "pr_1week")
    mslp = get(simdir, "mslp_1week")
    ts = get(simdir, "ts_1week")

    pr.data .*= 7
    pr.attributes["long_name"] = "Precipitation, accumulation within 1 Week"

    
    preprocessed_vars = JLD2.load_object("experiments/calibration/era5_preprocessed_vars.jld2")

    era5_pr = preprocessed_vars[findfirst(v -> ClimaAnalysis.short_name(v) == "pr", preprocessed_vars)]
    era5_ts = preprocessed_vars[findfirst(v -> ClimaAnalysis.short_name(v) == "tas", preprocessed_vars)]
    era5_mslp = preprocessed_vars[findfirst(v -> ClimaAnalysis.short_name(v) == "mslp", preprocessed_vars)]

    era5_mslp.attributes["long_name"] = "Mean sea level pressure, average within 1 Week"
    era5_ts.attributes["long_name"] = "Surface temperature, average within 1 Week"

    var_pairs = (
        (pr, era5_pr),
        (mslp, era5_mslp),
        (ts, era5_ts)
    )
    fig = GeoMakie.Figure(size = (1000, 500 * length(var_pairs)))
    for (i, (var, era5_var)) in enumerate(var_pairs)
        last_date = last(dates(var))
        var = slice(var, time = last_date)
        era5_var = slice(era5_var, time = last_date)
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[i,1], var - era5_var; more_kwargs = Dict(:plot => ca_kwargs(colormap = :redsblues)))
    end
    GeoMakie.save(joinpath(simdir.simulation_path, "bias_redsblues.png"), fig)
end