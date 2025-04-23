using ClimaAnalysis
using Statistics, Dates, LinearAlgebra, OrderedCollections, NaNStatistics

function auto_covariance(output_var)
    lon, lat = size(output_var.data)[1:2]
    covariance = zeros(Float32, lon, lat)
    for i in 1:lon
        for j in 1:lat
            covariance[i, j] = var(output_var.data[i, j, :])
        end
    end
    return covariance
end

function compute_covariance(var; full = true)
    monthly_cre = split_by_months(var)
    if !full # hard coding months here
        monthly_cre = [monthly_cre[9]]
    end
    monthly_covariances = vec.(auto_covariance.(monthly_cre))
    covariance_vec = vcat(monthly_covariances...)
    size_of_vec = length(covariance_vec)
    covariance_mat = Diagonal(covariance_vec)

    days_in_secs = 86_400
    monthly_time_avgs = [var.data for var in ClimaAnalysis.average_time.(monthly_cre)]
    discrepancy  = Diagonal(Float32.(0.05 * vcat(vec.(monthly_time_avgs)...))).^2

    # TODO: It might be helpful to make a function in ClimaAnalysis to vectorize this
    # cre_2010 = window(var, "time"; left = 93days_in_secs, right = 426days_in_secs)
    # discrepancy = Diagonal(Float32.(vec(cre_2010.data)))

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

# Move the two functions below to ClimaAnalysis
function split_by_months(var::OutputVar)
    ClimaAnalysis.Var._check_time_dim(var)
    start_date = Dates.DateTime(var.attributes["start_date"])

    monthly_dates = split_by_months(ClimaAnalysis.Utils.time_to_date.(start_date, ClimaAnalysis.Var.times(var)))
    monthly_times =
        (ClimaAnalysis.Var.date_to_time.(start_date, month) for month in monthly_dates)

    monthly_vars = ClimaAnalysis.Var._split_along_dim(var, ClimaAnalysis.Var.time_name(var), monthly_times)
    return monthly_vars
end

function split_by_months(dates::AbstractArray{<:Dates.DateTime})
    # Dates are not necessarily sorted
    dates = sort(dates)

    # Empty case
    isempty(dates) && return Vector{Vector{eltype(dates)}}[]

    month2dates = OrderedDict{Int, Vector{eltype(dates)}}()
    for i in 1:12
        month2dates[i] = []
    end

    for date in dates
        month = Dates.month(date)
        push!(month2dates[month], date)
    end

    return collect(values(month2dates))
end
