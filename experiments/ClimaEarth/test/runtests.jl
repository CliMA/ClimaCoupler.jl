using SafeTestsets
import ClimaComms
ClimaComms.@import_required_backends

gpu_broken = ClimaComms.device() isa ClimaComms.CUDADevice

@safetestset "component model test: ClimaAtmos" begin
    include("component_model_tests/climaatmos_tests.jl")
end
@safetestset "component model test: ClimaLand integrated model" begin
    include("component_model_tests/climaland_tests.jl")
end
@safetestset "component model test: prescr. sea ice" begin
    include("component_model_tests/prescr_seaice_tests.jl")
end
@safetestset "component model test: prescr. sea ice" begin
    include("component_model_tests/prescr_ocean_tests.jl")
end
gpu_broken || @safetestset "component model test: eisenman sea ice" begin
    include("component_model_tests/eisenman_seaice_tests.jl")
end
@safetestset "component model test: slab ocean" begin
    include("component_model_tests/slab_ocean_tests.jl")
end
@safetestset "debug diagnostics: debug plots" begin
    include("debug_plots_tests.jl")
end
@safetestset "AMIP test" begin
    include("amip_test.jl")
end

@safetestset "Restart test" begin
    include("restart.jl")
end
