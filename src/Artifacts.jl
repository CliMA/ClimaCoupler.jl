module Artifacts

using ClimaUtilities
using LazyArtifacts

"""
    precipitation_obs_data(; context = nothing)

Return the path for NetCDF file containing the observed monthly averaged precipitation data.
"""
function precipitation_obs_data(; context = nothing)
    return joinpath(
        ClimaUtilities.ClimaArtifacts.@clima_artifact("precipitation_obs", context),
        "gpcp.precip.mon.mean.197901-202305.nc",
    )
end

end
