const OBS_DS = Dict()
const SIM_DS_KWARGS = Dict()
const OTHER_MODELS_RMSEs = Dict()
const COMPARISON_RMSEs = Dict()

function preprocess_pr_fn(data)
    # -1 kg/m/s2 -> 1 mm/day
    return data .* Float32(-86400)
end

Base.@kwdef struct RMSEs
    model_name::String
    ANN::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
    DJF::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
    JJA::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
    MAM::Union{<:Real, ClimaAnalysis.OutputVar} = 0.0
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

include("cmip_rmse.jl")

# OBS_DS["rsut"] = ObsDataSource(;
#                              path = "OBS/CERES_EBAF-TOA_Ed4.2_Subset_200003-202303.g025.nc",
#                              var_name = "toa_sw_all_mon",
#                              )

# OBS_DS["rlut"] = ObsDataSource(;
#                              path = "OBS/CERES_EBAF-TOA_Ed4.2_Subset_200003-202303.g025.nc",
#                              var_name = "toa_lw_all_mon",
#                              )

function bias(output_dir::AbstractString, short_name::AbstractString, target_dates::AbstractArray{<:Dates.DateTime})
    obs = OBS_DS[short_name]
    sim = SimDataSource(; path = output_dir, short_name, SIM_DS_KWARGS["pr"]...)
    return bias(obs, sim, target_dates)
end

function compute_biases(output_dir, short_names, target_dates::AbstractArray{<:Dates.DateTime})
    return map(name -> bias(output_dir, name, target_dates), short_names)
end

function plot_biases(biases; output_path)
    fig = CairoMakie.Figure(; size = (600, 300 * length(biases)))
    loc = 1
    for bias_var in biases
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig, bias_var; p_loc = (1, loc))
        loc = loc + 1
    end
    CairoMakie.save(output_path, fig)
end

function plot_leaderboard(rmses; output_path)
    fig = CairoMakie.Figure(; size = (600, 300 * length(rmses)))
    loc = 1

    for rmse in rmses
        short_name = rmse.ANN.attributes["var_short_name"]
        units = rmse.ANN.attributes["units"]
        ax = CairoMakie.Axis(
            fig[1, loc],
            ylabel = "$short_name [$units]",
            xticks = (1:5, ["Ann", "DJF", "JJA", "MAM", "SON"]),
            title = "Global RMSE",
        )

        # Against other models
        (; best_single_model, median_model, worst_model, best_model) = COMPARISON_RMSEs[short_name]

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
    CairoMakie.save(output_path, fig)
end
