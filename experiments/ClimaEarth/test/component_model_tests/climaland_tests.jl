using Test
import ClimaCore as CC
import ClimaAtmos as CA
import ClimaCoupler
import ClimaCoupler: FluxCalculator, Interfacer
import Dates
import Thermodynamics.Parameters as TDP
import ClimaParams # to load TDP extension
import ClimaComms
ClimaComms.@import_required_backends

exp_dir = joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth")

include(joinpath(exp_dir, "components", "land", "climaland_integrated.jl"))
include(joinpath(exp_dir, "components", "atmosphere", "climaatmos.jl"))

FT = Float32

@testset "ClimaLandSimulation constructor" begin
    dt = Float64(450)
    tspan = (Float64(0), 3.0dt)
    start_date = Dates.DateTime(2008)
    output_dir = pwd()
    boundary_space = CC.CommonSpaces.CubedSphereSpace(
        FT;
        radius = FT(6371e3),
        n_quad_points = 4,
        h_elem = 4,
    )
    area_fraction = CC.Fields.ones(boundary_space)
    atmos_h = CC.Fields.zeros(boundary_space) .+ 2

    # Construct simulation object
    land_sim =
        ClimaLandSimulation(FT; dt, tspan, start_date, output_dir, area_fraction, atmos_h)

    # Try taking a timestep
    Interfacer.step!(land_sim, dt)

    # Check that the simulation object is correctly initialized
    @test nameof(land_sim) == "ClimaLandSimulation"
    @test land_sim.area_fraction == area_fraction

    # Check that the state is correctly initialized
    state_names = propertynames(land_sim.integrator.u)
    @test :canopy in state_names
    @test :soil in state_names
    @test :snow in state_names
    @test :soilco2 in state_names

    # Check that the cache is correctly initialized
    cache_names = propertynames(land_sim.integrator.p)
    @test :canopy in cache_names
    @test :soil in cache_names
    @test :snow in cache_names
    @test :soilco2 in cache_names
    @test :drivers in cache_names

    # Check that the drivers are correctly initialized
    driver_names = propertynames(land_sim.integrator.p.drivers)
    @test driver_names ==
          (:P_liq, :P_snow, :c_co2, :T, :P, :q, :SW_d, :LW_d, :cosθs, :frac_diff, :soc)
end

@testset "ClimaLandSimulation flux calculations" begin
    dt = Float64(120)
    tspan = (Float64(0), 3.0dt)
    start_date = Dates.DateTime(2008)
    output_dir = pwd()

    # Construct atmos and land simulation objects
    atmos_config_file =
        joinpath(exp_dir, "test", "component_model_tests", "climaatmos_coarse_short.yml")
    atmos_config = CA.AtmosConfig(atmos_config_file; job_id = "atmos_land_flux_test")
    atmos_sim = ClimaAtmosSimulation(atmos_config)

    boundary_space = ClimaCore.Spaces.horizontal_space(atmos_sim.domain.face_space)
    area_fraction = ClimaCore.Fields.ones(boundary_space)
    atmos_h = CC.Fields.zeros(boundary_space) .+ 2
    land_sim =
        ClimaLandSimulation(FT; dt, tspan, start_date, output_dir, area_fraction, atmos_h)
    model_sims = (; land_sim = land_sim, atmos_sim = atmos_sim)

    # Initialize the coupler fields so we can perform exchange
    coupler_field_names = Interfacer.default_coupler_fields()
    map(sim -> Interfacer.add_coupler_fields!(coupler_field_names, sim), values(model_sims))
    coupler_fields = Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)
    thermo_params = TDP.ThermodynamicsParameters(FT)

    cs = Interfacer.CoupledSimulation{FT}(
        nothing, # start_date
        coupler_fields,
        nothing, # conservation_checks
        tspan,
        dt,
        tspan[1],
        Ref(-1), # prev_checkpoint_t
        model_sims,
        (;), # callbacks
        (;), # dir_paths
        thermo_params, # thermo_params
        nothing, # diags_handler
    )

    # Step the atmosphere once to get non-zero wind and humidity so the fluxes are non-zero
    Interfacer.step!(atmos_sim, dt)

    # Exchange the initial conditions between atmosphere and land
    # This also tests the `get_field`, `update_field!` and `update_model_sims!` methods for `ClimaLandSimulation`
    FieldExchanger.exchange!(cs)

    # Update land cache variables with the updated drivers in the cache after the exchange
    update_aux! = CL.make_update_aux(land_sim.model)
    update_aux!(land_sim.integrator.p, land_sim.integrator.u, land_sim.integrator.t)

    update_boundary_fluxes! = CL.make_update_boundary_fluxes(land_sim.model)
    update_boundary_fluxes!(
        land_sim.integrator.p,
        land_sim.integrator.u,
        land_sim.integrator.t,
    )

    # Compute the surface fluxes
    FluxCalculator.compute_surface_fluxes!(
        coupler_fields,
        land_sim,
        atmos_sim,
        thermo_params,
    )

    # Check that the fluxes have been changed
    zero_field = CC.Fields.zeros(boundary_space)
    @test coupler_fields.F_turb_ρτxz != zero_field
    @test coupler_fields.F_turb_ρτyz != zero_field
    @test coupler_fields.F_lh != zero_field
    @test coupler_fields.F_sh != zero_field
    @test coupler_fields.F_turb_moisture != zero_field

    # Check that the fluxes don't contain any NaNs
    @test !any(isnan, coupler_fields.F_turb_ρτxz)
    @test !any(isnan, coupler_fields.F_turb_ρτyz)
    @test !any(isnan, coupler_fields.F_lh)
    @test !any(isnan, coupler_fields.F_sh)
    @test !any(isnan, coupler_fields.F_turb_moisture)

    # Check that drivers in cache got updated
    for driver in propertynames(land_sim.integrator.p.drivers)
        # Snow and liquid precipitation are zero with this setup
        if !(driver in [:P_liq, :P_snow])
            @test getproperty(land_sim.integrator.p.drivers, driver) != zero_field
        end
    end
end
