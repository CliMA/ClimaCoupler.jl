using ClimaAnalysis
import Dates
import ClimaCalibrate
import GeoMakie
import Makie
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import ClimaCoupler
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

    for m in 1:EKP.get_N_ens(ekp)
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
    if !ClimaCalibrate.EnsembleBuilder.is_complete(g_ens_builder)
        for m in 1:EKP.get_N_ens(ekp)
            mean(g_ens_builder.g_ens[:, m]) ≈ 0 && (g_ens_builder.g_ens[:, m] .= NaN)
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
        short_name = short_name * "_1week" # TODO: Don't hardcode this
        var = get(simdir, short_name)
        var.attributes["short_name"] = replace(short_name, "_1week" => "")

        var = preprocess_var(var, sample_date_ranges[iteration+1] )

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
function preprocess_var(var::ClimaAnalysis.OutputVar, sample_date_range)
    # TODO: Check for sign of pr
    # TODO: Check for units of everything
    var = shift_to_previous_week(var)
    @assert ClimaAnalysis.short_name(var) in CALIBRATE_CONFIG.short_names

    if ClimaAnalysis.short_name(var) in ("pr", "tas")
        var = ClimaAnalysis.apply_oceanmask(var)
    end
    
    var = window(var, "time"; left = sample_date_range[1], right = sample_date_range[2])
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
    simdir = SimDir(ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1))

    plot_bias(simdir; output_dir = plot_output_path)
    plot_variables(simdir; output_dir = plot_output_path)
    plot_spread_in_variables(iteration, ekp)
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

function plot_variables(simdir; output_dir = simdir.simulation_path)
    vars = map(x -> get(simdir, x * "_1week"), CALIBRATE_CONFIG.short_names)
    fig = GeoMakie.Figure(size = (1000, length(vars)*500))

    for (i, var) in enumerate(vars)
        var = slice(var, time = last(dates(var)))
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[i, 1], var; more_kwargs = Dict(:plot => ca_kwargs(colormap = :viridis)))
    end

    GeoMakie.save(joinpath(output_dir, "vars.png"), fig)
end

function plot_bias(simdir; output_dir = simdir.simulation_path)
    vars = Dict()
    for short_name in CALIBRATE_CONFIG.short_names
        vars[short_name] = get(simdir, short_name * "_1week")
    end
    # TODO: Get observations from EKP object
    preprocessed_vars = JLD2.load_object("experiments/calibration/era5_preprocessed_vars.jld2")

    era5_pr = preprocessed_vars[findfirst(v -> ClimaAnalysis.short_name(v) == "pr", preprocessed_vars)]
    era5_tas = preprocessed_vars[findfirst(v -> ClimaAnalysis.short_name(v) == "tas", preprocessed_vars)]
    era5_mslp = preprocessed_vars[findfirst(v -> ClimaAnalysis.short_name(v) == "mslp", preprocessed_vars)]

    era5_mslp.attributes["long_name"] = "Mean sea level pressure, average within 1 Week"
    era5_tas.attributes["long_name"] = "Surface temperature, average within 1 Week"

    var_pairs = (
        (vars["pr"], era5_pr),
        (vars["mslp"], era5_mslp),
        (vars["tas"], era5_tas)
    )
    plot_extrema = Dict(
        "tas" => (-3, 3), 
        "mslp" => (-1000, 1000),
        "pr" => (-1e-4, 1e-4),
    )

    fig = GeoMakie.Figure(size = (1000, 500 * length(var_pairs)))
    for (i, (var, era5_var)) in enumerate(var_pairs)
        last_date = last(dates(var))
        var = slice(var, time = last_date)
        era5_var = slice(era5_var, time = last_date)
        global_bias = ClimaAnalysis.global_bias(var, era5_var)
        global_mean = weighted_average_lonlat(var).data[1]
        @info " $(short_name(var)): Global bias = $global_bias, mean = $global_mean"

        global_bias = ClimaAnalysis.global_bias(var, era5_var; mask = ClimaAnalysis.apply_oceanmask)
        land_mean = weighted_average_lonlat(ClimaAnalysis.apply_oceanmask(var)).data[1]
        @info " $(short_name(var)): Land bias = $global_bias, mean = $land_mean"

        global_bias = ClimaAnalysis.global_bias(var, era5_var; mask = ClimaAnalysis.apply_landmask)
        ocean_mean = weighted_average_lonlat(ClimaAnalysis.apply_landmask(var)).data[1]
        @info " $(short_name(var)): Ocean bias = $global_bias, mean = $ocean_mean"
        ClimaAnalysis.Visualize.plot_bias_on_globe!(fig[i,1], var, era5_var; cmap_extrema = plot_extrema[short_name(var)])
    end
    GeoMakie.save(joinpath("bias_redsblues.png"), fig)
end

function absolute_global_bias(sim_var, obs_var)
    global_bias = abs(ClimaAnalysis.global_bias(sim_var, obs_var).data[1])
    global_mean = abs(ClimaAnalysis.average_lonlat(sim_var).data[1])
    return (global_bias, global_mean)
end

function plot_spread_in_variables(iteration, ekp)
    short_names = CALIBRATE_CONFIG.short_names
    var_spread = Dict(short_name => [] for short_name in short_names)
    for m in 1:EKP.get_N_ens(ekp)
        member_path = ClimaCalibrate.path_to_ensemble_member(CALIBRATE_CONFIG.output_dir, iteration, m)
        simdir = SimDir(member_path)
        for short_name in short_names
            try
                var = get(simdir, short_name * "_1week")
                var_spread[short_name] = [var_spread[short_name]..., var]
            catch e
                @warn e
            end
        end
    end
    fig = GeoMakie.Figure(size = (1000, length(short_names) * 500))
    for (i, short_name) in enumerate(short_names)
        ensemble_data = getproperty.(var_spread[short_name], :data)
        plot_var = remake(var_spread[short_name][1]; data = std(ensemble_data))
        last_slice = slice(plot_var, time = last(dates(plot_var)))
        # TODO: Don't just use the last slice
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[i,1], last_slice; more_kwargs = Dict( :axis => ca_kwargs(title = "Standard Deviation $(short_name)")))
        # TODO: Set titles
    end
    GeoMakie.save(joinpath(ClimaCalibrate.path_to_iteration(CALIBRATE_CONFIG.output_dir, iteration), "ensemble_stdev.png"), fig)

end