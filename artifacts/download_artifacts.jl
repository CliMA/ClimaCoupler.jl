
# Download atmos artifacts
module Atmos
import ClimaAtmos
include(joinpath(pkgdir(ClimaAtmos), "artifacts", "download_artifacts.jl"))
trigger_download()
end

# Download land artifacts
module Land
import ClimaLand
if pkgversion(ClimaLand) < v"0.15.2"
    include(joinpath(pkgdir(ClimaLand), "src", "standalone", "Bucket", "artifacts", "artifacts.jl"))
    @info "CESM", cesm2_albedo_dataset_path()
    @info "Bareground albedo", bareground_albedo_dataset_path()
end
end


include(joinpath(@__DIR__, "artifact_funcs.jl"))

# Trigger download if data doesn't exist locally
function trigger_download()
    @info "co2 dataset path: `$(co2_dataset_path())`"
    @info "mask dataset path: `$(mask_dataset_path())`"
    return nothing
end
trigger_download()
