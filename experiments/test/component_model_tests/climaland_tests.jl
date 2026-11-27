using Test
import ClimaCore as CC
import ClimaAtmos as CA
import ClimaCoupler
import ClimaCoupler: FieldExchanger, FluxCalculator, Interfacer
import Dates
import Thermodynamics.Parameters as TDP
import ClimaParams as CP # to load TDP extension
import ClimaComms
ClimaComms.@import_required_backends

exp_dir = joinpath(pkgdir(ClimaCoupler), "experiments")

import ClimaAtmos
# Needed to construct ClimaAtmosSimulation
ClimaAtmosExt = Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt)

# To load ClimaCouplerClimaLandExt
import ClimaLand as CL
import NCDatasets

FT = Float32

@testset "ClimaLandSimulation constructor" begin
    coupled_param_dict = CP.create_toml_dict(FT)

    dt = Float64(450)
    tspan = (Float64(0), 3.0dt)
    start_date = Dates.DateTime(2008)
    output_dir = pwd()
    boundary_space = CC.CommonSpaces.CubedSphereSpace(
        FT;
        radius = coupled_param_dict["planet_radius"], # in meters
        n_quad_points = 4,
        h_elem = 4,
    )
    area_fraction = CC.Fields.ones(boundary_space)
    atmos_h = CC.Fields.zeros(boundary_space) .+ 2
    coord = CC.Fields.coordinate_field(boundary_space)
    initial_T = CC.Fields.zeros(boundary_space) .+ 276.85 .+ 40 .* cosd.(coord.lat) .^ 4
    # Construct simulation object
    land_sim = Interfacer.LandSimulation(
        FT,
        Val(:integrated);
        dt,
        tspan,
        start_date,
        output_dir,
        area_fraction,
        surface_space = boundary_space,
        atmos_h,
        initial_T,
    )

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
          (:P_liq, :P_snow, :c_co2, :T, :P, :q, :u, :SW_d, :LW_d, :cosθs, :frac_diff)
    atmos = land_sim.model.soil.boundary_conditions.top.atmos
    @test atmos == land_sim.model.canopy.boundary_conditions.atmos
    @test atmos == land_sim.model.snow.boundary_conditions.atmos
    # Remap atmos_h to the same space as atmos.h for type comparison
    atmos_h_remapped = Interfacer.remap(axes(atmos.h), atmos_h)
    @test atmos.h == atmos_h_remapped
    #@test typeof(atmos.h) == typeof(atmos_h_remapped)
    @test atmos.gustiness == FT(1)
    # Newer ClimaLand versions also store the height `reference` convention
    @test issubset((:h, :gustiness), propertynames(atmos))
end

@testset "ClimaLandSimulation flux calculations" begin
    coupled_param_dict = CP.create_toml_dict(FT)
    thermo_params = TDP.ThermodynamicsParameters(coupled_param_dict)

    dt = Float64(120)
    tspan = (Float64(0), 3.0dt)
    start_date = Dates.DateTime(2008)
    output_dir = pwd()

    # Construct atmos and land simulation objects
    atmos_config_file =
        joinpath(exp_dir, "test", "component_model_tests", "climaatmos_coarse_short.yml")
    atmos_config = CA.AtmosConfig(atmos_config_file; job_id = "atmos_land_flux_test")
    atmos_sim = ClimaAtmosExt.ClimaAtmosSimulation(atmos_config)

    boundary_space = CC.Spaces.horizontal_space(atmos_sim.domain.face_space)
    area_fraction = CC.Fields.ones(boundary_space)
    atmos_h = CC.Fields.zeros(boundary_space) .+ 2
    initial_T = Interfacer.get_field(boundary_space, atmos_sim, Val(:air_temperature))
    land_sim = Interfacer.LandSimulation(
        FT,
        Val(:integrated);
        dt,
        tspan,
        start_date,
        output_dir,
        area_fraction,
        surface_space = boundary_space,
        atmos_h,
        initial_T,
    )
    model_sims = (; land_sim = land_sim, atmos_sim = atmos_sim)

    # Initialize the coupler fields so we can perform exchange
    coupler_field_names = Interfacer.default_coupler_fields()
    map(sim -> Interfacer.add_coupler_fields!(coupler_field_names, sim), values(model_sims))
    coupler_fields = Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)

    cs = Interfacer.CoupledSimulation{FT}(
        nothing, # start_date
        coupler_fields,
        nothing, # conservation_checks
        tspan,
        dt,
        Ref(tspan[1]),
        Ref(0), # step
        Ref(-1), # prev_checkpoint_t
        model_sims,
        (;), # callbacks
        (;), # dir_paths
        thermo_params, # thermo_params
        nothing, # diags_handler
        true, # save_cache
        Interfacer.SteppingSchedule(false, false), # schedule
        (;), # flux_accumulators
    )

    # Step the atmosphere once to get non-zero wind and humidity so the fluxes are non-zero
    Interfacer.step!(atmos_sim, dt)

    # Exchange the initial conditions between atmosphere and land
    # This also tests the `get_field`, `update_field!` and `update_model_sims!` methods for `ClimaLandSimulation`
    FieldExchanger.import_static_fields!(coupler_fields, model_sims)
    FieldExchanger.exchange!(cs)

    # Initialize land cache after drivers are exchanged. This sets
    # `p.snow.T_sfc` (used as the SurfaceFluxes ≥1.2 initial guess for the snow
    # skin-temperature Newton update) and updates aux / boundary-flux cache
    # variables.
    Interfacer.set_cache!(land_sim, coupler_fields)

    # Compute the turbulent fluxes for each sub-component as the land model does,
    # with the atmospheric height set by `set_cache!`; the soil and snow fluxes
    # come from their skin temperature solves, not from `CL.turbulent_fluxes!`
    update_boundary_fluxes! = CL.make_update_boundary_fluxes(land_sim.model)
    update_boundary_fluxes!(
        land_sim.integrator.p,
        land_sim.integrator.u,
        land_sim.integrator.t,
    )

    # Combine the surface fluxes from each sub-component and update the coupler fields
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

@testset "ClimaLand long-period diagnostic compute subsampling" begin
    import ClimaDiagnostics as CD
    ClimaLandExt = Base.get_extension(ClimaCoupler, :ClimaCouplerClimaLandExt)

    coupled_param_dict = CP.create_toml_dict(FT)
    start_date = Dates.DateTime(2008, 1, 1)
    boundary_space = CC.CommonSpaces.CubedSphereSpace(
        FT;
        radius = coupled_param_dict["planet_radius"],
        n_quad_points = 4,
        h_elem = 2,
    )
    area_fraction = CC.Fields.ones(boundary_space)
    atmos_h = CC.Fields.zeros(boundary_space) .+ FT(2)
    initial_T = CC.Fields.zeros(boundary_space) .+ FT(280)

    function scheduled_diags(sim)
        for cb in sim.integrator.callback.discrete_callbacks
            hasproperty(cb.affect!, :diagnostics_handler) &&
                return cb.affect!.diagnostics_handler.scheduled_diagnostics
        end
        error("No diagnostics_handler found")
    end

    mktempdir() do outdir
        # Monthly average with dt = 450 s (< 3600 s) subsamples hourly, and a
        # mid-hour restart (`tspan[1] = 900 s`) keeps calendar-hour phase.
        bucket_monthly = Interfacer.LandSimulation(
            FT,
            Val(:bucket);
            dt = 450.0,
            tspan = (900.0, 4500.0),
            start_date,
            output_dir = joinpath(outdir, "monthly"),
            area_fraction,
            surface_space = boundary_space,
            atmos_h,
            initial_T,
            albedo_type = "function",
            land_diagnostics_period = :monthly,
            land_diagnostics_reduction = :average,
        )
        sds_monthly = scheduled_diags(bucket_monthly)
        @test !isempty(sds_monthly)
        @test all(sds_monthly) do sd
            sd.compute_schedule_func isa CD.Schedules.EveryCalendarDtSchedule &&
                sd.compute_schedule_func.dt == Dates.Hour(1) &&
                sd.compute_schedule_func.date_last[] == start_date
        end
        # Stateful schedules are distinct per diagnostic.
        @test sds_monthly[1].compute_schedule_func !== sds_monthly[2].compute_schedule_func

        # Coarse dt >= 3600 s keeps EveryStepSchedule.
        bucket_coarse_dt = Interfacer.LandSimulation(
            FT,
            Val(:bucket);
            dt = 3600.0,
            tspan = (0.0, 7200.0),
            start_date,
            output_dir = joinpath(outdir, "coarse_dt"),
            area_fraction,
            surface_space = boundary_space,
            atmos_h,
            initial_T,
            albedo_type = "function",
            land_diagnostics_period = :monthly,
            land_diagnostics_reduction = :average,
        )
        @test all(scheduled_diags(bucket_coarse_dt)) do sd
            sd.compute_schedule_func isa CD.Schedules.DivisorSchedule &&
                sd.compute_schedule_func.divisor == 1
        end

        # Daily average keeps EveryStepSchedule.
        bucket_daily = Interfacer.LandSimulation(
            FT,
            Val(:bucket);
            dt = 450.0,
            tspan = (0.0, 1800.0),
            start_date,
            output_dir = joinpath(outdir, "daily"),
            area_fraction,
            surface_space = boundary_space,
            atmos_h,
            initial_T,
            albedo_type = "function",
            land_diagnostics_period = :daily,
            land_diagnostics_reduction = :average,
        )
        @test all(scheduled_diags(bucket_daily)) do sd
            sd.compute_schedule_func isa CD.Schedules.DivisorSchedule &&
                sd.compute_schedule_func.divisor == 1
        end
    end
end
