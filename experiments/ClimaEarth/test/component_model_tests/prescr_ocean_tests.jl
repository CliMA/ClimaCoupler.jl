using Test
import Dates
import ClimaCore as CC
import Thermodynamics.Parameters as TDP
import ClimaParams # required for TDP
import ClimaCoupler

include(joinpath("..", "TestHelper.jl"))
import .TestHelper
include(joinpath("..", "..", "components", "ocean", "prescr_ocean.jl"))

FT = Float32

@testset "PrescribedOceanSimulation name" begin
    sim = PrescribedOceanSimulation((;))
    @test Interfacer.name(sim) == "PrescribedOceanSimulation"
end

@testset "PrescribedOceanSimulation constructor" begin
    space = TestHelper.create_space(FT)
    date0 = Dates.DateTime(2000, 1, 1)
    t_start = 0.0
    area_fraction = CC.Fields.ones(space)
    thermo_params = TDP.ThermodynamicsParameters(FT)
    comms_ctx = nothing

    # Construct simulation object
    sim = PrescribedOceanSimulation(FT, space, date0, t_start, area_fraction, thermo_params, comms_ctx)

    # Read in initial SST data
    sst_data = try
        joinpath(@clima_artifact("historical_sst_sic", comms_ctx), "MODEL.SST.HAD187001-198110.OI198111-202206.nc")
    catch error
        @warn "Using lowres SST. If you want the higher resolution version, you have to obtain it from ClimaArtifacts"
        joinpath(
            @clima_artifact("historical_sst_sic_lowres", comms_ctx),
            "MODEL.SST.HAD187001-198110.OI198111-202206_lowres.nc",
        )
    end

    SST_timevaryinginput = TimeVaryingInput(
        sst_data,
        "SST",
        space,
        reference_date = date0,
        file_reader_kwargs = (; preprocess_func = (data) -> data + FT(273.15),), ## convert to Kelvin
    )
    SST_expected = zeros(space)
    evaluate!(SST_expected, SST_timevaryinginput, t_start)

    @test sim isa Interfacer.AbstractSurfaceStub
    # Check that the cache is correctly initialized
    @test sim.cache.T_sfc == SST_expected
    @test sim.cache.ρ_sfc == zeros(space)
    @test sim.cache.z0m == FT(5.8e-5)
    @test sim.cache.z0b == FT(5.8e-5)
    @test sim.cache.beta == FT(1)
    @test sim.cache.α_direct == ones(space) .* FT(0.06)
    @test sim.cache.α_diffuse == ones(space) .* FT(0.06)
    @test sim.cache.area_fraction == area_fraction
    @test sim.cache.phase == TD.Liquid()
    @test sim.cache.thermo_params == thermo_params
    @test !isnothing(sim.cache.SST_timevaryinginput)

    # Test `Interfacer.get_field` function
    @test Interfacer.get_field(sim, Val(:air_density)) == sim.cache.ρ_sfc
    @test Interfacer.get_field(sim, Val(:area_fraction)) == sim.cache.area_fraction
    @test Interfacer.get_field(sim, Val(:beta)) == sim.cache.beta
    @test Interfacer.get_field(sim, Val(:roughness_buoyancy)) == sim.cache.z0b
    @test Interfacer.get_field(sim, Val(:roughness_momentum)) == sim.cache.z0m
    @test Interfacer.get_field(sim, Val(:surface_direct_albedo)) == sim.cache.α_direct
    @test Interfacer.get_field(sim, Val(:surface_diffuse_albedo)) == sim.cache.α_diffuse
    @test Interfacer.get_field(sim, Val(:surface_temperature)) == sim.cache.T_sfc
end
