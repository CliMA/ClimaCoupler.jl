import SafeTestsets: @safetestset
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
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
gpu_broken || @safetestset "Regridder tests" begin
    include("regridder_tests.jl")
end
@safetestset "ConservationChecker tests" begin
    include("conservation_checker_tests.jl")
end
gpu_broken || @safetestset "BCReader tests" begin
    include("bcreader_tests.jl")
end
@safetestset "Utilities tests" begin
    include("utilities_tests.jl")
end
@safetestset "TimeManager tests" begin
    include("time_manager_tests.jl")
end
@safetestset "FieldExchanger tests" begin
    include("field_exchanger_tests.jl")
end
@safetestset "FluxCalculator tests" begin
    include("flux_calculator_tests.jl")
end
gpu_broken || @safetestset "Diagnostics tests" begin
    include("diagnostics_tests.jl")
end
gpu_broken || @safetestset "PostProcessor tests" begin
    include("postprocessor_tests.jl")
end
@safetestset "Checkpointer tests" begin
    include("checkpointer_tests.jl")
end
gpu_broken || @safetestset "experiment test: CoupledSims tests" begin
    include("experiment_tests/coupled_sims.jl")
end
@safetestset "experiment test: Leaderboard" begin
    include("experiment_tests/leaderboard.jl")
end
@safetestset "component test: bucket" begin
    include("component_model_tests/bucket_tests.jl")
end
@safetestset "component model test: ClimaAtmos" begin
    include("component_model_tests/climaatmos_tests.jl")
end
@safetestset "component model test: prescr. sea ice" begin
    include("component_model_tests/prescr_seaice_tests.jl")
end
gpu_broken || @safetestset "component model test: eisenman sea ice" begin
    include("component_model_tests/eisenman_seaice_tests.jl")
end
@safetestset "component model test: slab ocean" begin
    include("component_model_tests/slab_ocean_tests.jl")
end
@safetestset "debug diagnostics: amip plots" begin
    include("debug/debug_amip_plots.jl")
end
