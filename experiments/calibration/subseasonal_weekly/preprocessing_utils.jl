import ClimaAnalysis

"""
    apply_lat_window(vars, lat_left, lat_right)

Apply latitude window by constraining the longitudes to be in the range
[lat_left, lat_right].
"""
function apply_lat_window(vars::AbstractDict, lat_left, lat_right)
    result = empty(vars)
    for (key, var) in vars
        lats = ClimaAnalysis.latitudes(var)
        first_lat_idx = findfirst(l -> l >= lat_left, lats)
        last_lat_idx = findlast(l -> l <= lat_right, lats)
        var = ClimaAnalysis.window(var, "latitude", by = ClimaAnalysis.Index(), left = first_lat_idx, right = last_lat_idx)
        @info "Latitudes of $(ClimaAnalysis.short_name(var)) is $(ClimaAnalysis.latitudes(var))"
        result[key] = var
    end
    return result
end

function apply_lat_window(vars::AbstractVector, lat_left, lat_right)
    preprocessed_vars = ClimaAnalysis.OutputVar[]
    for var in vars
        lats = ClimaAnalysis.latitudes(var)
        first_lat_idx = findfirst(l -> l >= lat_left, lats)
        last_lat_idx = findlast(l -> l <= lat_right, lats)
        var = ClimaAnalysis.window(var, "latitude", by = ClimaAnalysis.Index(), left = first_lat_idx, right = last_lat_idx)
        @info "Latitudes of $(ClimaAnalysis.short_name(var)) is $(ClimaAnalysis.latitudes(var))"
        push!(preprocessed_vars, var)
    end
    return preprocessed_vars
end
