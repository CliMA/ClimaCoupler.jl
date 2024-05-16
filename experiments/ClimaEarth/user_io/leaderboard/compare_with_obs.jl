const OBS_DS = Dict()

function preprocess_pr_fn(data)
    # 1 mm/day -> - 1 kg/m/s2
    # The minus sign comes from the different conventions used
    return data .* Float32(-1 / 86400)
end

OBS_DS["pr"] = ObsDataSource(;
    path = joinpath(pr_obs_data_path(), "gpcp.precip.mon.mean.197901-202305.nc"),
    var_name = "precip",
    preprocess_data_fn = preprocess_pr_fn,
)

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
    sim = SimDataSource(; path = output_dir, short_name)
    return bias(obs, sim, target_dates)
end

function plot_biases(output_dir, short_names, target_dates::AbstractArray{<:Dates.DateTime}; output_path)
    fig = CairoMakie.Figure(; size = (600, 300 * length(short_names)))
    loc = 1
    for short_name in short_names
        bias_var = bias(output_dir, short_name, target_dates)
        ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig, bias_var; p_loc = (1, loc))
        loc = loc + 1
    end
    CairoMakie.save(output_path, fig)
end
