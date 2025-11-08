using Statistics
using ClimaAnalysis
import Dates
import ClimaCalibrate
import GeoMakie
import Makie
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import ClimaCoupler
import ClimaCalibrate: EnsembleBuilder
using OrderedCollections

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
            EnsembleBuilder.fill_g_ens_col!(g_ens_builder, m, NaN)
        end
    end
    if count(isnan, g_ens_builder.g_ens) > 0.9 * length(g_ens_builder.g_ens)
        error("Too many NaNs")
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
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges[iteration+1]
    @info "Short names: $short_names"

    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in short_names
        var = get_var(short_name, simdir)
        var = preprocess_var(var, sample_date_ranges )
        @show dates(var)

        EnsembleBuilder.fill_g_ens_col!(g_ens_builder, col_idx, var; verbose = true)
    end

    return nothing
end

const period = "1M"

function get_var(short_name, simdir)
    latent_heat_vaporization_at_reference = 2500800  # J kg^-1

    # TODO: detect period from the sample date range
    if short_name == "tas - ta"
        tas = get(simdir; short_name = "tas")
        ta = get(simdir; short_name = "ta")
        # pfull = get(simdir; short_name = "pfull")
        # ta_900 = ClimaAnalysis.Atmos.to_pressure_coordinates(ta, pfull; target_pressure=[900])
        ta_900 = slice(ta; z = 1000)
        var = tas - ta_900
    elseif short_name == "tas"
        var = get(simdir; short_name = "tas")
    elseif short_name == "rlns"
        rlus = get(simdir; short_name = "rlus")
        rlds = get(simdir; short_name = "rlds")
        var = rlus - rlds
    elseif short_name == "rsns"
        rsus = get(simdir; short_name = "rsus")
        rsds = get(simdir; short_name = "rsds")
        var = rsus - rsds
    elseif short_name == "hfls"
        var = convert_units(get(simdir, "evspsbl"), "W m^-2"; conversion_function = x -> x * latent_heat_vaporization_at_reference)
        var.attributes["long_name"] = "Latent Heat Flux at the surface, average within 1 month"
    elseif short_name == "hfss"
        lhf = convert_units(get(simdir, "evspsbl"), "W m^-2"; conversion_function = x -> x * latent_heat_vaporization_at_reference)
        lhf.attributes["long_name"] = "Latent Heat Flux at the surface, average within 1 month"
        var = get(simdir, "hfes") - lhf
        var.attributes["long_name"] = "Sensible Heat Flux at the surface, average within 1 month"
    else
    
        # TODO: support multiple periods/reductions of same variable
        var = get(simdir; short_name)
    end

    var.attributes["short_name"] = short_name
    return var
end

"""
    preprocess_var(var::ClimaAnalysis.OutputVar, reference_date)

Preprocess `var` before flattening for G ensemble matrix.

For "pr", weekly sums are computed. For "tas" and "mslp", weekly means are
computed from daily means. The daily means are computing starting from
`reference_date`.

This function assumes that the data is monthly.
"""
function preprocess_var(var, sample_date_range)
    @show short_name(var)
    @assert short_name(var) in CALIBRATE_CONFIG.short_names
    var = if period == "1M"
       shift_to_start_of_previous_month(var)
    elseif period == "1W"
       shift_to_previous_week(var)
    end
    if short_name(var) in ("sw_cre","hfss", "hfls", "rlns", "rsns")
        var = set_units(var, "W m^-2")
    elseif short_name(var) in ("tas", "tas - ta", "ta")
        var = set_units(var, "K")
    end
    @show units(var)

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

    simdir = SimDir(ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, 1))
    plot_bias(ekp, simdir, iteration; output_dir = plot_output_path)
    # plot_variables(simdir; output_dir = plot_output_path)
    # plot_pointwise_spread_per_variable(ekp, iteration)
    try 
        plot_surface_fluxes_bias(simdir; output_dir = plot_output_path)
    catch e
        @error e
    end

    @info "Ensemble spread: $(scalar_spread(ekp, iteration))" 
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

"""
    plot_bias(ekp, simdir, iteration; output_dir)

Plot bias maps comparing simulation output to observations for all variables in
`CALIBRATE_CONFIG.short_names`.

Uses observations from the EKP object via `reconstruct_vars` and compares
them against simulation variables for each date in the sample date ranges.
"""
function plot_bias(ekp, simdir, iteration; output_dir = simdir.simulation_path)
    # Get observations for this iteration
    obs_series = EKP.get_observation_series(ekp)
    minibatch_obs = ClimaCalibrate.ObservationRecipe.get_observations_for_nth_iteration(
        obs_series,
        iteration + 1
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
        var = shift_to_start_of_previous_month(var)
        # Apply unit conversions
        if short_name in ("tas", "tas - ta", "ta")
            var = set_units(var, "K")
        elseif short_name in ("hfls", "hfss", "rsns", "rlns")
            var = set_units(var, "W m^-2")
        end
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
    
    # Color map extrema for different variables
    plot_extrema = Dict(
        "tas" => (-6, 6), 
        "tas - ta" => (-6, 6), 
        "hfls" => (-50, 50),
        "hfss" => (-25, 25),
        "rsns" => (-50, 50),
        "rlns" => (-50, 50),
        "mslp" => (-1000, 1000),
        "pr" => (-1e-4, 1e-4),
    )
    
    # Get sample dates for this iteration
    sample_dates = unique(CALIBRATE_CONFIG.sample_date_ranges[iteration + 1])
    
    # Create figure
    fig = GeoMakie.Figure(size = (1500, 500 * length(var_pairs)))
    
    for (i, (sim_var, era5_var)) in enumerate(var_pairs)
        for (j, date) in enumerate(sample_dates)
            sim_var_t = slice(sim_var, time = date)
            era5_var_t = slice(era5_var, time = date)
            
            # Calculate biases
            global_bias = ClimaAnalysis.global_bias(sim_var_t, era5_var_t)
            global_mean = weighted_average_lonlat(sim_var_t).data[1]
            relative_global_bias = global_bias / global_mean
            
            land_bias = ClimaAnalysis.global_bias(sim_var_t, era5_var_t; mask = ClimaAnalysis.apply_oceanmask)
            land_mean = weighted_average_lonlat(ClimaAnalysis.apply_oceanmask(sim_var_t)).data[1]
            relative_land_bias = land_bias / land_mean
            
            ocean_bias = ClimaAnalysis.global_bias(sim_var_t, era5_var_t; mask = ClimaAnalysis.apply_landmask)
            ocean_mean = weighted_average_lonlat(ClimaAnalysis.apply_landmask(sim_var_t)).data[1]
            relative_ocean_bias = ocean_bias / ocean_mean
            
            @info short_name(sim_var_t) relative_global_bias global_bias global_mean
            @info short_name(sim_var_t) relative_land_bias land_bias land_mean
            @info short_name(sim_var_t) relative_ocean_bias ocean_bias ocean_mean
            
            # Get color map extrema
            cmap_extrema = get(plot_extrema, short_name(sim_var_t), extrema(sim_var_t.data))
            @show cmap_extrema
            
            # Plot bias
            ClimaAnalysis.Visualize.plot_bias_on_globe!(fig[i, j], sim_var_t, era5_var_t; cmap_extrema)
        end
    end
    
    GeoMakie.save(joinpath(output_dir, "bias_sample_dates.png"), fig)
    return nothing
end



function plot_pointwise_spread_per_variable(ekp, iteration)
    ensemble_vars = get_ensemble_of_vars(ekp, iteration)
    fig = GeoMakie.Figure(size = (1000, length(CALIBRATE_CONFIG.short_names) * 500))
    for (i, short_name) in enumerate(CALIBRATE_CONFIG.short_names)
        ensemble_data = getproperty.(ensemble_vars[short_name], :data)
        data = std(ensemble_data; corrected = false)
        plot_var = remake(ensemble_vars[short_name][1]; data)
        last_slice = slice(plot_var, time = last(dates(plot_var)))
        # TODO: Don't just use the last slice
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[i,1], last_slice; more_kwargs = Dict( :axis => ca_kwargs(title = "Stdev $(short_name)")))
        # TODO: Set titles
    end
    GeoMakie.save(joinpath(ClimaCalibrate.path_to_iteration(CALIBRATE_CONFIG.output_dir, iteration), "ensemble_stdev.png"), fig)

end


function scalar_spread(ekp, iteration)
    g_mean_final = EKP.get_g_mean_final(ekp)
    g_final = EKP.get_g_final(ekp)
    sq_dists = [sum((col .- g_mean_final).^2) for col in eachcol(g_final)]
    return mean(sq_dists) 
end

function get_ensemble_of_vars(ekp, iteration)
    short_names = CALIBRATE_CONFIG.short_names
    ensemble_vars = Dict(short_name => [] for short_name in short_names)
    for m in 1:EKP.get_N_ens(ekp)
        member_path = ClimaCalibrate.path_to_ensemble_member(CALIBRATE_CONFIG.output_dir, iteration, m)
        simdir = SimDir(member_path)
        for short_name in short_names
            try
                var = get(simdir, short_name * "_1week")
                ensemble_vars[short_name] = [ensemble_vars[short_name]..., var]
            catch e
                @warn e
            end
        end
    end
    return ensemble_vars
end


cmap_extrema = Dict(
    "lhf" => (0, 220),
    "shf" => (0, 80 ),
    "rsds - rsus" => (0, 300),
    "rlds - rlus" => (-120, 0)
)

function plot_surface_fluxes(simdir; output_dir = simdir.simulation_path)

    latent_heat_vaporization_at_reference = 2500800  # J kg^-1
    # LHF units: W m^-2 = ( J  kg^-1 ) * ( kg m^-2 s^-1 )
    lhf = convert_units(get(simdir, "evspsbl"), "W m^-2"; conversion_function = x -> x * latent_heat_vaporization_at_reference)
    lhf.attributes["long_name"] = "Latent Heat Flux at the surface, average within 1 month"
    lhf.attributes["short_name"] = "lhf"
    shf = get(simdir, "hfes") - lhf
    shf.attributes["long_name"] = "Sensible Heat Flux at the surface, average within 1 month"
    shf.attributes["short_name"] = "shf"

    rsds = get(simdir, "rsds")
    rsus = get(simdir, "rsus")
    
    rsns = rsds - rsus
    rsns.attributes["long_name"] = "Net SW Flux at the surface, average within 1 month"

    rlds = get(simdir, "rlds")
    rlus = get(simdir, "rlus")
    rlns = rlds - rlus
    rlns.attributes["long_name"] = "Net LW Flux at the surface, average within 1 month"

    vars = [lhf, shf, rsns, rlns]

    fig = GeoMakie.Figure(size = (1000, length(vars)*500))
    for (i, var) in enumerate(vars)
        var = average_time(var)
        min_val, max_val = cmap_extrema[short_name(var)]
        
        # Create levels between min and max
        levels = collect(range(min_val, max_val, length = 11))  # 11 levels
        
        ClimaAnalysis.Visualize.contour2D_on_globe!(
            fig[i, 1], 
            var;
            more_kwargs = Dict(
                :plot => merge(
                    ca_kwargs(colormap = :viridis, extendhigh = :auto, extendlow = :auto),
                    Dict(:levels => levels)
                )
            )
        )
    end

    GeoMakie.save(joinpath(output_dir, "surface_fluxes.png"), fig)

end

function preprocess_era5(var_name)
    short_name_map = Dict(
        "avg_slhtf" => "lhf",
        "avg_ishf" => "shf",
        "avg_snswrf" => "rsns",
        "avg_snlwrf" => "rlns",
    )
    var = OutputVar("era5_monthly_mean_surface_fluxes.nc", var_name)
    reverse_dim!(var, latitude_name(var))
    var.attributes["short_name"] = short_name_map[var_name]
    if short_name(var) == "lhf" || short_name(var) == "shf"
        var = -(var)
        var.attributes["short_name"] = short_name_map[var_name]

    end
    var = set_units(var, "W m^-2")
    # Get rid of missing type
    data = ones(size(var.data)...)
    data .= var.data
    dims = OrderedDict([dim => Float64[var.dims[dim]...] for dim in keys(var.dims)])
    var = remake(var; data, dims)
    return var
end


function plot_surface_fluxes_bias(simdir; output_dir = joinpath(simdir.simulation_path,"../.."))

    bias_cmap_extrema = Dict(
        "lhf" => (-50, 50),
        "shf" => (-25, 25, ),
        "rsds - rsus" => (-50, 50),
        "rlds - rlus" => (-50, 50)
    )
    latent_heat_vaporization_at_reference = 2500800  # J kg^-1
    # LHF units: W m^-2 = ( J  kg^-1 ) * ( kg m^-2 s^-1 )
    lhf = convert_units(get(simdir, "evspsbl"), "W m^-2"; conversion_function = x -> x * latent_heat_vaporization_at_reference)
    lhf.attributes["long_name"] = "Latent Heat Flux at the surface, average within 1 month"
    lhf.attributes["short_name"] = "lhf"
    shf = get(simdir, "hfes") - lhf
    shf.attributes["long_name"] = "Sensible Heat Flux at the surface, average within 1 month"
    shf.attributes["short_name"] = "shf"
    shf = set_units(shf, "W m^-2")

    rsds = get(simdir, "rsds")
    rsus = get(simdir, "rsus")
    
    rsns = rsds - rsus
    rsns.attributes["long_name"] = "Net SW Flux at the surface, average within 1 month"
    rsns = set_units(rsns, "W m^-2")

    rlds = get(simdir, "rlds")
    rlus = get(simdir, "rlus")
    rlns = rlds - rlus
    rlns.attributes["long_name"] = "Net LW Flux at the surface, average within 1 month"
    rlns = set_units(rlns, "W m^-2")
    
    vars_sim = [lhf, shf, rsns, rlns]

    vars_sim = shift_to_start_of_previous_month.(vars_sim)

    vars_era5 = preprocess_era5.(["avg_slhtf", "avg_ishf", "avg_snswrf", "avg_snlwrf"])
    start_date = DateTime(2018,9,1)

    vars_sim = window.(vars_sim, "time", left = start_date, right = start_date)
    vars_era5 = window.(vars_era5, "time", left = start_date, right = start_date)
    @info dates(first(vars_era5))
    @info dates(first(vars_sim))
    fig = GeoMakie.Figure(size = (1000, length(vars_sim)*500));
    for (i, (var_sim, var_era5)) in enumerate(zip(vars_sim, vars_era5))

        var_sim = average_time(var_sim)
        var_era5 = average_time(var_era5)
        @show short_name(var_era5)
        min_val, max_val = bias_cmap_extrema[short_name(var_sim)]
        
        # Create levels between min and max
        levels = collect(range(min_val, max_val, length = 11))  # 11 levels
        
        ClimaAnalysis.Visualize.plot_bias_on_globe!(
            fig[i, 1], 
            var_sim,
            var_era5;
            cmap_extrema = (min_val, max_val),
            more_kwargs = Dict(
                :plot => merge(
                    ca_kwargs(colormap = :vik, extendhigh = :auto, extendlow = :auto),
                    Dict(:levels => levels)
                )
            )
        )
    end

    GeoMakie.save(joinpath(output_dir, "surface_fluxes_bias.png"), fig)

end

