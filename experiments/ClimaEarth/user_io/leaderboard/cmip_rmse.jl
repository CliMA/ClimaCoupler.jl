import Statistics: median, quantile

const RMSE_FILE_PATHS = Dict()

RMSE_FILE_PATHS["pr"] = joinpath(@clima_artifact("cmip_model_rmse"), "pr_rmse_amip_pr_amip_5yr.csv")
RMSE_FILE_PATHS["rsut"] = joinpath(@clima_artifact("cmip_model_rmse"), "rsut_rmse_amip_rsut_amip_5yr.csv")
RMSE_FILE_PATHS["rlut"] = joinpath(@clima_artifact("cmip_model_rmse"), "rlut_rmse_amip_rlut_amip_5yr.csv")

short_names = ["pr", "rsut", "rlut"]
for short_name in short_names
    open(RMSE_FILE_PATHS[short_name], "r") do io
        # Skip the header
        header = readline(io)

        # Process each line
        for line in eachline(io)
            # Split the line by comma
            fields = split(line, ',')
            model_name = fields[1]
            DJF, MAM, JJA, SON, ANN = map(x -> parse(Float64, x), fields[2:end])

            push!(OTHER_MODELS_RMSEs[short_name], RMSEs(; model_name, DJF, MAM, JJA, SON, ANN))
        end
    end
end

"""
    best_single_model(RMSEs)

Return the one model that has the overall smallest error.
"""
function best_single_model(RMSEs)
    _, index = findmin(r -> abs.(values(r)), RMSEs)
    return RMSEs[index]
end

"""
    worst_single_model(RMSEs)

Return the one model that has the overall largest error.
"""
function worst_single_model(RMSEs)
    _, index = findmax(r -> abs.(values(r)), RMSEs)
    return RMSEs[index]
end

"""
    RMSE_stats(RMSEs)

RMSEs is the dictionary OTHER_MODELS_RMSEs.

Return:
- best single model
- worst single model
- "model" with all the medians
- "model" with all the best values
- "model" with all the worst values
"""
function RMSE_stats(dict_vecRMSEs)
    stats = Dict()
    # cumulative_error maps model_names with the total RMSE across metrics normalized by median(RMSE)
    cumulative_error = Dict()
    for (key, vecRMSEs) in dict_vecRMSEs
        # Collect into vectors that we can process independently
        all_values = stack(values.(vecRMSEs))
        ANN, DJF, JJA, MAM, SON = ntuple(i -> all_values[i, :], 5)

        median_model = RMSEs(;
            model_name = "Median",
            ANN = median(ANN),
            DJF = median(DJF),
            JJA = median(JJA),
            MAM = median(MAM),
            SON = median(SON),
        )

        worst_model = RMSEs(;
            model_name = "Worst",
            ANN = maximum(abs.(ANN)),
            DJF = maximum(abs.(DJF)),
            JJA = maximum(abs.(JJA)),
            MAM = maximum(abs.(MAM)),
            SON = maximum(abs.(SON)),
        )

        best_model = RMSEs(;
            model_name = "Best",
            ANN = minimum(abs.(ANN)),
            DJF = minimum(abs.(DJF)),
            JJA = minimum(abs.(JJA)),
            MAM = minimum(abs.(MAM)),
            SON = minimum(abs.(SON)),
        )

        quantile25 = RMSEs(;
            model_name = "Quantile 0.25",
            ANN = quantile(ANN, 0.25),
            DJF = quantile(DJF, 0.25),
            JJA = quantile(JJA, 0.25),
            MAM = quantile(MAM, 0.25),
            SON = quantile(SON, 0.25),
        )

        quantile75 = RMSEs(;
            model_name = "Quantile 0.75",
            ANN = quantile(ANN, 0.75),
            DJF = quantile(DJF, 0.75),
            JJA = quantile(JJA, 0.75),
            MAM = quantile(MAM, 0.75),
            SON = quantile(SON, 0.75),
        )

        for rmse in vecRMSEs
            haskey(cumulative_error, cumulative_error) || (cumulative_error[rmse.model_name] = 0.0)
            cumulative_error[rmse.model_name] += sum(values(rmse) ./ values(median_model))
        end

        stats[key] = (;
            best_single_model = best_single_model(vecRMSEs),
            worst_single_model = worst_single_model(vecRMSEs),
            median_model,
            worst_model,
            best_model,
            quantile25,
            quantile75,
        )
    end

    _, absolute_best_model = findmin(cumulative_error)
    _, absolute_worst_model = findmax(cumulative_error)

    return (; stats, absolute_best_model, absolute_worst_model)
end

const COMPARISON_RMSEs_STATS = RMSE_stats(OTHER_MODELS_RMSEs)
