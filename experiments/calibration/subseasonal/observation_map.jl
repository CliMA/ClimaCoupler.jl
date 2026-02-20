using Statistics
import JLD2
import Dates
using ClimaAnalysis
import ClimaCalibrate
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import ClimaCoupler
import ClimaCalibrate: EnsembleBuilder

include(joinpath(@__DIR__, "observation_utils.jl"))

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
            @error "Ensemble member $m failed" exception = (e, catch_backtrace())
            EnsembleBuilder.fill_g_ens_col!(g_ens_builder, m, NaN)
        end
    end
    g_ens = EnsembleBuilder.get_g_ensemble(g_ens_builder)
    if count(isnan, g_ens) > 0.9 * length(g_ens)
        error("Too many NaNs")
    end
    return EnsembleBuilder.is_complete(g_ens_builder) ? g_ens :
           error("G ensemble matrix is not completed")
end

"""
    process_member_data!(diagnostics_folder_path, short_names, current_minibatch)

Process the data of a single ensemble member and return a single column of the
G ensemble matrix.
"""
function process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx, iteration)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges[iteration + 1]
    @info "Short names: $short_names"

    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in short_names
        var = get_var(short_name, simdir)
        var = preprocess_var(var, sample_date_ranges)

        EnsembleBuilder.fill_g_ens_col!(
            g_ens_builder,
            col_idx,
            var;
            checkers = (EnsembleBuilder.SequentialIndicesChecker(),),
            verbose = true,
        )
    end

    return nothing
end

function largest_period(sample_date_range)
    span = maximum(sample_date_range) - minimum(sample_date_range)
    span = Millisecond(span)
    span.value == 0 && return Month(1)
    day_in_ms = 8.64e7
    period =
        span.value >= day_in_ms * 365 ? Year(1) :
        span.value >= day_in_ms * 30 ? Month(1) :
        span.value >= day_in_ms * 7 ? Week(1) : Day(1)
    return period
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

    simdir = SimDir(ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1))
    plot_bias(ekp, simdir, iteration; output_dir = plot_output_path)

    @info "Ensemble spread: $(scalar_spread(ekp))"
end

"""
    plot_constrained_params_and_errors(output_dir, ekp, prior)

Plot the constrained parameters and errors from `ekp` and `prior` and save
them to `output_dir`.
"""
function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 1], ekp, error_metric = "loss")
    EKP.Visualize.plot_error_over_time(fig[1, dim_size + 2], ekp, error_metric = "loss")
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
end

bias_plot_extrema = Dict(
    "tas" => (-6, 6),
    "tas - ta" => (-6, 6),
    "hfls" => (-50, 50),
    "hfss" => (-25, 25),
    "rsus" => (-50, 50),
    "rlus" => (-50, 50),
    "mslp" => (-1000, 1000),
    "pr" => (-1e-4, 1e-4),
)

"""
    plot_bias(ekp, simdir, iteration; output_dir)

Plot bias maps comparing simulation output to observations for all variables in
`CALIBRATE_CONFIG.short_names`.

Uses observations from the EKP object via `reconstruct_vars` and compares
them against simulation variables for each date in the sample date ranges.
"""
function plot_bias(ekp, simdir, iteration; output_dir = simdir.simulation_path)
    # Get observations for this iteration
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges[iteration + 1]
    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs = ClimaCalibrate.ObservationRecipe.get_observations_for_nth_iteration(
        obs_series,
        iteration + 1,
    )

    # Reconstruct OutputVars from observations (ERA5 data)
    era5_vars = []
    for obs in minibatch_obs
        obs_vars = ClimaCalibrate.ObservationRecipe.reconstruct_vars(obs)
        append!(era5_vars, obs_vars)
    end

    # Get simulation variables for all short_names
    sim_vars = []
    for short_name in CALIBRATE_CONFIG.short_names
        var = get_var(short_name, simdir)
        var = preprocess_var(var, sample_date_ranges)
        push!(sim_vars, var)
    end

    # Match sim_vars with era5_vars by short_name
    var_pairs = []
    for sim_var in sim_vars
        sim_short_name = ClimaAnalysis.short_name(sim_var)
        era5_idx = findfirst(v -> ClimaAnalysis.short_name(v) == sim_short_name, era5_vars)
        if !isnothing(era5_idx)
            push!(var_pairs, (sim_var, era5_vars[era5_idx]))
        else
            @warn "No ERA5 data found for $(sim_short_name)"
        end
    end

    if isempty(var_pairs)
        @warn "No matching variable pairs found for bias plotting"
        return nothing
    end

    # Get sample dates for this iteration
    sample_dates = unique(CALIBRATE_CONFIG.sample_date_ranges[iteration + 1])

    # Create figure
    fig = GeoMakie.Figure(size = (1500, 500 * length(var_pairs)))
    for (j, date) in enumerate(sample_dates)
        for (i, (sim_var, era5_var)) in enumerate(var_pairs)
            sim_var_t = slice(sim_var, time = date)
            era5_var_t = slice(era5_var, time = date)

            # Calculate biases
            global_bias = ClimaAnalysis.global_bias(sim_var_t, era5_var_t)
            global_mean = weighted_average_lonlat(sim_var_t).data[1]
            relative_global_bias = global_bias / global_mean

            land_bias = ClimaAnalysis.global_bias(
                sim_var_t,
                era5_var_t;
                mask = ClimaAnalysis.apply_oceanmask,
            )
            land_mean =
                weighted_average_lonlat(ClimaAnalysis.apply_oceanmask(sim_var_t)).data[1]
            relative_land_bias = land_bias / land_mean

            ocean_bias = ClimaAnalysis.global_bias(
                sim_var_t,
                era5_var_t;
                mask = ClimaAnalysis.apply_landmask,
            )
            ocean_mean =
                weighted_average_lonlat(ClimaAnalysis.apply_landmask(sim_var_t)).data[1]
            relative_ocean_bias = ocean_bias / ocean_mean

            @info short_name(sim_var_t) relative_global_bias global_bias global_mean
            @info short_name(sim_var_t) relative_land_bias land_bias land_mean
            @info short_name(sim_var_t) relative_ocean_bias ocean_bias ocean_mean

            cmap_extrema =
                get(bias_plot_extrema, short_name(sim_var_t), extrema(sim_var_t.data))

            # Plot bias
            ax = ClimaAnalysis.Visualize.plot_bias_on_globe!(
                fig[i, j],
                sim_var_t,
                era5_var_t;
                cmap_extrema,
            )
        end
    end

    GeoMakie.save(joinpath(output_dir, "bias_sample_dates.png"), fig)
    return nothing
end


function scalar_spread(ekp)
    g_mean_final = EKP.get_g_mean_final(ekp)
    g_final = EKP.get_g_final(ekp)
    sq_dists = [sum((col .- g_mean_final) .^ 2) for col in eachcol(g_final)]
    return mean(sq_dists)
end
