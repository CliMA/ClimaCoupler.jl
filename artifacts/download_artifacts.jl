
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
