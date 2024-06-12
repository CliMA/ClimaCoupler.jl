const OBS_DS = Dict()
const SIM_DS_KWARGS = Dict()
const OTHER_MODELS_RMSEs = Dict()
const COMPARISON_RMSEs = Dict()

function preprocess_pr_fn(data)
    # -1 kg/m/s2 -> 1 mm/day
    return data .* Float32(-86400)
end

function replace_nan(x; v = 0.0)
    return map(x -> isnan(x) ? zero(x) : x, v)
end

Base.@kwdef struct RMSEs
    model_name::String
    ANN::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
    DJF::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
    MAM::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
    JJA::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
    SON::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
end

function Base.values(r::RMSEs)
    val_or_rmse(v::Real) = v
    val_or_rmse(v::ClimaAnalysis.OutputVar) = v.attributes["rmse"]

    return val_or_rmse.([r.ANN, r.DJF, r.JJA, r.MAM, r.SON])
end

OBS_DS["pr"] = ObsDataSource(;
    path = joinpath(@clima_artifact("precipitation_obs"), "gpcp.precip.mon.mean.197901-202305.nc"),
    var_name = "precip",
)
SIM_DS_KWARGS["pr"] = (; preprocess_data_fn = preprocess_pr_fn, new_units = "mm / day")
OTHER_MODELS_RMSEs["pr"] = []

OBS_DS["rsut"] = ObsDataSource(;
    path = joinpath(@clima_artifact("radiation_obs"), "CERES_EBAF-TOA_Ed4.2_Subset_200003-202303.g025.nc"),
    var_name = "toa_sw_all_mon",
)
SIM_DS_KWARGS["rsut"] = (;)
OTHER_MODELS_RMSEs["rsut"] = []

OBS_DS["rlut"] = ObsDataSource(;
    path = joinpath(@clima_artifact("radiation_obs"), "CERES_EBAF-TOA_Ed4.2_Subset_200003-202303.g025.nc"),
    var_name = "toa_lw_all_mon",
)
SIM_DS_KWARGS["rlut"] = (;)
OTHER_MODELS_RMSEs["rlut"] = []

include("cmip_rmse.jl")

function bias(output_dir::AbstractString, short_name::AbstractString, target_dates::AbstractArray{<:Dates.DateTime})
    obs = OBS_DS[short_name]
    sim = SimDataSource(; path = output_dir, short_name, SIM_DS_KWARGS[short_name]...)
    return bias(obs, sim, target_dates)
end

function compute_biases(output_dir, short_names, target_dates::AbstractArray{<:Dates.DateTime})
    return map(name -> bias(output_dir, name, target_dates), short_names)
end

function plot_biases(biases; output_path)
    fig = CairoMakie.Figure(; size = (600, 300 * length(biases)))
    loc = 1
    for bias_var in biases
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig, bias_var; p_loc = (loc, 1))
        loc = loc + 1
    end
    CairoMakie.save(output_path, fig)
end

function plot_leaderboard(rmses; output_path)
    fig = CairoMakie.Figure(; size = (800, 300 * length(rmses) + 400), fontsize = 20)
    loc = 1

    NUM_BOXES = 4 + 1 # 4 seasons and 1 annual
    NUM_MODELS = 2  # CliMA vs best

    num_variables = length(rmses)

    var_names = map(r -> r.ANN.attributes["var_short_name"], rmses)

    # The square plot is at the very bottom
    loc_squares = length(rmses) + 1
    ax_squares = CairoMakie.Axis(
        fig[loc_squares, 1],
        yticks = (1:num_variables, var_names),
        xticks = ([3, NUM_BOXES + 3], ["CliMA", "Best model"]),
        aspect = NUM_BOXES * NUM_MODELS,
    )
    ax_squares2 = CairoMakie.Axis(
        fig[loc_squares, 1],
        xaxisposition = :top,
        xticks = (0.5:4.5, ["Ann", "DJF", "MAM", "JJA", "SON"]),
        aspect = NUM_BOXES * NUM_MODELS,
    )
    CairoMakie.hidespines!(ax_squares2)
    CairoMakie.hideydecorations!(ax_squares2)

    # Preallocate the matrix for the squares, each row is (4 seasons + annual) x
    # models compared, and there is one row per variable
    squares = zeros(NUM_BOXES * NUM_MODELS, num_variables)

    for (var_num, rmse) in enumerate(rmses)
        short_name = rmse.ANN.attributes["var_short_name"]
        units = rmse.ANN.attributes["units"]
        ax = CairoMakie.Axis(
            fig[loc, 1],
            ylabel = "$short_name [$units]",
            xticks = (1:5, ["Ann", "DJF", "MAM", "JJA", "SON"]),
            title = "Global RMSE $short_name [$units]",
        )

        # Against other models
        (; best_single_model, median_model, worst_model, best_model) = COMPARISON_RMSEs[short_name]

        squares[begin:NUM_BOXES, var_num] .= values(rmse) ./ values(median_model)
        squares[(NUM_BOXES + 1):end, var_num] .= values(best_single_model) ./ values(median_model)

        CairoMakie.errorbars!(
            ax,
            1:5,
            values(median_model),
            values(best_model),
            values(worst_model),
            whiskerwidth = 10,
            color = :black,
            linewidth = 0.5,
        )
        CairoMakie.scatter!(
            ax,
            1:5,
            values(median_model),
            label = median_model.model_name,
            color = :black,
            marker = :hline,
        )
        CairoMakie.scatter!(ax, 1:5, values(best_single_model), label = best_single_model.model_name)
        CairoMakie.scatter!(ax, 1:5, values(rmse), label = rmse.model_name, marker = :star5)

        # Add a fake extra point to center the legend a little better
        CairoMakie.scatter!(ax, [6.5], [0.1], markersize = 0.01)
        CairoMakie.axislegend()
        loc = loc + 1
    end

    colormap = CairoMakie.Reverse(:RdYlGn)

    # Now, the square plot
    CairoMakie.heatmap!(
        ax_squares,
        squares,
        colormap = colormap,
        # Trick to exclude the zeros
        lowclip = :white,
        colorrange = (1e-10, maximum(squares)),
    )
    CairoMakie.vlines!(ax_squares2, NUM_BOXES, color = :black, linewidth = 3.0)
    CairoMakie.Colorbar(
        fig[loc_squares, 2],
        limits = extrema(squares),
        label = "RMSE/median(RMSE)",
        colormap = colormap,
    )

    CairoMakie.save(output_path, fig)
end
