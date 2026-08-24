# This file contains tests for functions that can be used without the ClimaCouplerClimaAtmosExt
# extension, but have different behavior when the extension is not loaded because the
# extension defines additonal defaults.

using Test
import ClimaCoupler
import ClimaCoupler: Input, Interfacer, Models

if isnothing(Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt))
    @testset "atmos_default_config_dict and set_albedos! with and without ClimaAtmos loaded" begin
        @test_logs (:warn, r"ClimaAtmos") @test Input.atmos_default_config_dict() == Dict()
        # the error thrown for a PrescribedOceanSimulation should tell the user to load ClimaAtmos or define their own set_albedos!
        @test_logs (:error, r"ClimaAtmos") Interfacer.set_albedos!(
            Models.PrescribedOceanSimulation(nothing),
            nothing,
        )
        fallback_method =
            which(Interfacer.set_albedos!, (Models.PrescribedOceanSimulation, Nothing))
        # # this fallback should not tell the user to load ClimaAtmos, since it is not specific to PrescribedOceanSimulation
        @test_logs (:error, r"SurfaceStub") Interfacer.set_albedos!(
            Interfacer.SurfaceStub(nothing),
            nothing,
        )

        # # now load in ClimaAtmos and make sure the correct method is called instead of the fallback
        import ClimaAtmos
        atmos_default = Input.atmos_default_config_dict()
        @test atmos_default isa Dict && !isempty(atmos_default)
        # test that the method being called when ClimaAtmos is loaded is different
        atmos_default_set_albedos! =
            which(Interfacer.set_albedos!, (Models.PrescribedOceanSimulation, Nothing))
        @test atmos_default_set_albedos! != fallback_method
    end
else
    @info "ClimaCouplerClimaAtmosExt is loaded, skipping tests for fallbacks used when the extension is not loaded."
end
