import CairoMakie

const OBS_DS = Dict()
const SIM_DS_KWARGS = Dict()
const OTHER_MODELS_RMSEs = Dict()

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

    return val_or_rmse.([r.ANN, r.DJF, r.MAM, r.JJA, r.SON])
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

OBS_DS["rsutcs"] = ObsDataSource(;
    path = joinpath(@clima_artifact("radiation_obs"), "CERES_EBAF-TOA_Ed4.2_Subset_200003-202303.g025.nc"),
    var_name = "toa_sw_clr_c_mon",
)
SIM_DS_KWARGS["rsutcs"] = (;)

OBS_DS["rsdt"] = ObsDataSource(;
    path = joinpath(@clima_artifact("radiation_obs"), "CERES_EBAF-TOA_Ed4.2_Subset_200003-202303.g025.nc"),
    var_name = "solar_mon",
)
SIM_DS_KWARGS["rsdt"] = (;)

OBS_DS["rlutcs"] = ObsDataSource(;
    path = joinpath(@clima_artifact("radiation_obs"), "CERES_EBAF-TOA_Ed4.2_Subset_200003-202303.g025.nc"),
    var_name = "toa_lw_clr_c_mon",
)
SIM_DS_KWARGS["rlutcs"] = (;)

include("cmip_rmse.jl")

function bias(output_dir::AbstractString, short_name::AbstractString, target_dates::AbstractArray{<:Dates.DateTime})
    obs = OBS_DS[short_name]
    sim = SimDataSource(; path = output_dir, short_name, SIM_DS_KWARGS[short_name]...)
    return bias(obs, sim, target_dates)
end

function compute_biases(output_dir, short_names, target_dates::AbstractArray{<:Dates.DateTime}; cmap_extrema = Dict())
    return map(short_names) do name
        bias_outvar = bias(output_dir, name, target_dates)
        # The attribute is used in plot_bias to fix the colormap
        haskey(cmap_extrema, name) && (bias_outvar.attributes["cmap_extrema"] = cmap_extrema[name])
        return bias_outvar
    end
end


# colors are sampled in [0,1], so define a function that linearly transforms x ∈ [lo, hi] to [0, 1].
to_unitrange(x::Number, lo::Number, hi::Number) = (x - lo) / (hi - lo)

"""
    constrained_cmap(cols, lo, hi; [categorical=false], [rev=false], [mid=0])
    constrained_cmap(lo, hi; [categorical=false], [rev=false], [mid=0])

Constrain a colormap to a given range.

Given a colormap implicitly defined in `± maximum(abs, (lo, hi))`, constrain it to the range [lo, hi].
This is useful to ensure that a colormap which is desired to diverge symmetrically around zero maps
the same color intensity to the same magnitude.

The second form is a convenience function that uses the `:redsblues` colormap.

# Arguments
- `cols`: a vector of colors, or a ColorScheme
- `lo`: lower bound of the range
- `hi`: upper bound of the range

# Keyword Arguments
- `mid`: midpoint of the range  # TODO: test `mid` better
- `categorical`: flag for whether returned colormap should be categorical or continous
- `rev`: flag for whether to reverse the colormap

# Returns
- `cmap::Makie.ColorGradient`: a colormap
"""
function constrained_cmap(cols::Vector, lo, hi; mid = 0, categorical = false, rev = false)
    constrained_cmap(CairoMakie.Makie.ColorScheme(cols), lo, hi; mid, categorical, rev)
end
function constrained_cmap(cols::CairoMakie.Makie.ColorScheme, lo, hi; mid = 0, categorical = false, rev = false)
    rev && (cols = reverse(cols))  # reverse colorscheme if requested, don't reverse below in `cgrad`.
    absmax = maximum(abs, (lo, hi) .- mid)
    # map lo, hi ∈ [-absmax, absmax] onto [0,1] to sample their corresponding colors
    lo_m, hi_m = to_unitrange.((lo, hi) .- mid, -absmax, absmax)
    # values on [0,1] where each color in cols is defined
    colsvals = range(0, 1; length = length(cols))
    # filter colsvals, keep only values in [lo_m, hi_m] + the endpoints lo_m and hi_m.
    filter_colsvals = filter(x -> lo_m <= x <= hi_m, unique([lo_m; colsvals; hi_m]))
    # select colors in filtered range; interpolate new low and hi colors.
    newcols = CairoMakie.Makie.get(cols, filter_colsvals)
    # values on [0,1] where the new colors are defined
    new_colsvals = to_unitrange.(filter_colsvals, lo_m, hi_m)
    cmap = CairoMakie.cgrad(newcols, new_colsvals; categorical, rev = false)
    return cmap
end

function plot_biases(biases; output_path)
    fig = CairoMakie.Figure(; size = (600, 300 * length(biases)))
    loc = 1

    for bias_var in biases
        min_level, max_level = get(bias_var.attributes, "cmap_extrema", extrema(bias_var.data))

        # Make sure that 0 is at the center
        cmap = constrained_cmap(CairoMakie.cgrad(:vik).colors, min_level, max_level; categorical = true)
        nlevels = 11
        # Offset so that it covers 0
        levels = collect(range(min_level, max_level, length = nlevels))
        offset = levels[argmin(abs.(levels))]
        levels = levels .- offset
        ticklabels = map(x -> string(round(x; digits = 0)), levels)
        ticks = (levels, ticklabels)

        more_kwargs = Dict(
            :plot => Dict(:colormap => cmap, :levels => levels, :extendhigh => :auto, :extendlow => :auto),
            :cb => Dict(:ticks => ticks),
        )

        ClimaAnalysis.Visualize.contour2D_on_globe!(fig, bias_var; p_loc = (loc, 1), more_kwargs)
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
        yticks = (1:num_variables, reverse(var_names)),
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

    (; absolute_best_model, absolute_worst_model) = COMPARISON_RMSEs_STATS

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

        (; median_model) = COMPARISON_RMSEs_STATS.stats[short_name]

        best_single_model = first(filter(x -> x.model_name == absolute_best_model, OTHER_MODELS_RMSEs[short_name]))
        worst_single_model = first(filter(x -> x.model_name == absolute_worst_model, OTHER_MODELS_RMSEs[short_name]))

        squares[begin:NUM_BOXES, end - var_num + 1] .= values(rmse) ./ values(median_model)
        squares[(NUM_BOXES + 1):end, end - var_num + 1] .= values(best_single_model) ./ values(median_model)

        CairoMakie.scatter!(
            ax,
            1:5,
            values(median_model),
            label = median_model.model_name,
            color = :black,
            marker = :hline,
            markersize = 10,
            visible = false,
        )

        categories = vcat(map(_ -> collect(1:5), 1:length(OTHER_MODELS_RMSEs[short_name]))...)

        CairoMakie.boxplot!(
            ax,
            categories,
            vcat(values.(OTHER_MODELS_RMSEs[short_name])...);
            whiskerwidth = 1,
            width = 0.35,
            mediancolor = :black,
            color = :gray,
            whiskerlinewidth = 1,
        )

        CairoMakie.scatter!(ax, 1:5, values(best_single_model), label = absolute_best_model)
        CairoMakie.scatter!(ax, 1:5, values(worst_single_model), label = absolute_worst_model)

        # If we want to plot other models
        # for model in OTHER_MODELS_RMSEs[short_name]
        #     CairoMakie.scatter!(ax, 1:5, values(model), marker = :hline)
        # end

        CairoMakie.scatter!(
            ax,
            1:5,
            values(rmse),
            label = rmse.model_name,
            marker = :star5,
            markersize = 20,
            color = :green,
        )

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
