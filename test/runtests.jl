import SafeTestsets: @safetestset
import ClimaComms
ClimaComms.@import_required_backends
import Pkg, Artifacts

gpu_broken = ClimaComms.device() isa ClimaComms.CUDADevice

# Download test-only artifacts
#
# (Currently not natively supported by Julia)
artifacts_toml = joinpath(@__DIR__, "Artifacts.toml")
artifacts = Artifacts.select_downloadable_artifacts(artifacts_toml)
for name in keys(artifacts)
    Pkg.Artifacts.ensure_artifact_installed(name, artifacts[name], artifacts_toml)
end

@safetestset "Aqua tests" begin
    include("aqua.jl")
end
@safetestset "Interfacer tests" begin
    include("interfacer_tests.jl")
end

@safetestset "ConservationChecker tests" begin
    include("conservation_checker_tests.jl")
end
@safetestset "Utilities tests" begin
    include("utilities_tests.jl")
end
@safetestset "FieldExchanger tests" begin
    include("field_exchanger_tests.jl")
end
@safetestset "FluxCalculator tests" begin
    include("flux_calculator_tests.jl")
end
