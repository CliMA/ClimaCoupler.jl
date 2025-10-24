using ClimaAnalysis
import Dates
import ClimaCalibrate
import GeoMakie
import Makie
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import ClimaCoupler
import ClimaCalibrate: EnsembleBuilder


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
        simdir_path = joinpath(member_path, "amip_config/output_active")
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

    rsut = get(simdir, "rsut")
    rsutcs = get(simdir, "rsutcs")
    # var = rsutcs - rsut
    # var.attributes["short_name"] = "sw_cre"
    # var = set_units(var, "W m^-2")
    var = rsut
    var = shift_to_start_of_previous_month(var)
    dates = sample_date_ranges[iteration+1]
    # TODO: window var
    var = window(var, "time", left = dates[1], right = dates[2])
    EnsembleBuilder.fill_g_ens_col!(g_ens_builder, col_idx, var; verbose = true )

    return nothing
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
    plot_bias(simdir, iteration; output_dir = plot_output_path)
    plot_variables(simdir; output_dir = plot_output_path)
    plot_pointwise_spread_per_variable(ekp, iteration)
    # try 
    #     plot_surface_fluxes(simdir; output_dir = plot_output_path)
    # catch e
    #     @error e
    # end

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
    # vars = map(x -> get(simdir, x * "_1week"), CALIBRATE_CONFIG.short_names)
    # fig = GeoMakie.Figure(size = (1000, length(vars)*500))
    fig = GeoMakie.Figure(size = (1000, 500))
    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)

    rsut = get(simdir, "rsut")
    rsutcs = get(simdir, "rsutcs")
    # var = rsutcs - rsut
    # var.attributes["short_name"] = "sw_cre"
    var = rsut
    ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[i, 1], var; more_kwargs = Dict(:plot => ca_kwargs(colormap = :viridis)))
    GeoMakie.save(joinpath(output_dir, "vars.png"), fig)
end

# TODO: Generalize this 
function plot_bias(simdir, iteration; output_dir = simdir.simulation_path)
    var = get(simdir, "rsut") #- get(simdir, "rsut")
    var = set_units(var, "W m^-2")
    var = shift_to_start_of_previous_month(var)
    obs_var = JLD2.load_object(joinpath(pkgdir(ClimaCoupler),"experiments/calibration/era5_preprocessed_vars.jld2"))[1]
    fig = GeoMakie.Figure(size = (1500, 500))
    for (j, date) in enumerate(CALIBRATE_CONFIG.sample_date_ranges[iteration+1])
        var_t = slice(var, time = date)
        era5_var_t = slice(obs_var, time = date)
        global_bias = ClimaAnalysis.global_bias(var_t, era5_var_t)
        global_mean = weighted_average_lonlat(var_t).data[1]
        relative_global_bias = global_bias / global_mean

        land_bias = ClimaAnalysis.global_bias(var_t, era5_var_t; mask = ClimaAnalysis.apply_oceanmask)
        land_mean = weighted_average_lonlat(ClimaAnalysis.apply_oceanmask(var_t)).data[1]
        relative_land_bias = land_bias / land_mean

        ocean_bias = ClimaAnalysis.global_bias(var_t, era5_var_t; mask = ClimaAnalysis.apply_landmask)
        ocean_mean = weighted_average_lonlat(ClimaAnalysis.apply_landmask(var_t)).data[1]
        relative_ocean_bias = ocean_bias / ocean_mean
        @info short_name(var_t) relative_global_bias global_bias global_mean
        @info short_name(var_t) relative_land_bias land_bias land_mean
        @info short_name(var_t) relative_ocean_bias ocean_bias ocean_mean
        # relative_global_bias relative_land_bias relative_ocean_bias
        ClimaAnalysis.Visualize.plot_bias_on_globe!(fig[1, j], var_t, era5_var_t;)
    end
    GeoMakie.save(joinpath(ClimaCalibrate.path_to_iteration(CALIBRATE_CONFIG.output_dir, iteration), "bias_sample_dates.png"), fig)
end

using Statistics


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
    lhf.attributes["long_name"] = "Latent Heat Flux at the surface, avg 1mon"
    lhf.attributes["short_name"] = "lhf"
    shf = get(simdir, "hfes") - lhf
    shf.attributes["long_name"] = "Sensible Heat Flux at the surface, avg 1mon"
    shf.attributes["short_name"] = "shf"

    rsds = get(simdir, "rsds")
    rsus = get(simdir, "rsus")

    rsns = rsds - rsus
    rsns.attributes["long_name"] = "Net SW Flux at the surface, avg 1mon"

    rlds = get(simdir, "rlds")
    rlus = get(simdir, "rlus")
    rlns = rlds - rlus
    rlns.attributes["long_name"] = "Net LW Flux at the surface, avg 1mon"

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