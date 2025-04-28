using ClimaAnalysis
using Statistics, Dates, LinearAlgebra, OrderedCollections, NaNStatistics

function compute_covariance(var; full = true)
    monthly_cre = ClimaAnalysis.split_by_month(var)
    if !full # hard coding months here
        monthly_cre = [monthly_cre[9]]
    end
    monthly_covariances = vec.(Float32.(covariance_var.data) for covariance)
    covariance_vec = vcat(monthly_covariances...)
    size_of_vec = length(covariance_vec)
    covariance_mat = Diagonal(covariance_vec)

    monthly_time_avgs = [var.data for var in ClimaAnalysis.average_time.(monthly_cre)]
    discrepancy  = Diagonal(Float32.(0.05 * vcat(vec.(monthly_time_avgs)...))).^2

    regularization = (0.1)^2 .* Float32.(Diagonal(I, size_of_vec))

    lat_weights = Float32.(Diagonal(repeat(vec(lat_weight_mat(var)), length(monthly_cre))))
    return lat_weights .* (covariance_mat .+ discrepancy .+ regularization)
end

# Maybe move this function to ClimaAnalysis
function lat_weight_mat(output_var::OutputVar)
    lon_size, lat_size = size(output_var.data)[1:2]
    lat_weights = zeros(Float32, lon_size, lat_size)
    for i in 1:lon_size
        for j in 1:lat_size
            lat_weights[i, j] = 1.0f0 ./ max(cosd(ClimaAnalysis.latitudes(output_var)[j]), 0.1)
        end
    end
    return lat_weights
end
