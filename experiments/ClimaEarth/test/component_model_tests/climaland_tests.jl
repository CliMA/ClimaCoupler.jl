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

exp_dir = joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth")

import ClimaAtmos
# Needed to construct ClimaAtmosSimulation
ClimaAtmosExt = Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt)

include(joinpath(exp_dir, "setup_run.jl"))

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

    # Construct simulation object
    land_sim = Interfacer.LandSimulation(
        FT,
        Val(:integrated);
        dt,
        tspan,
        start_date,
        output_dir,
        area_fraction,
        atmos_h,
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
    @test propertynames(atmos) == (:h, :gustiness)
end

@testset "ClimaLandSimulation flux calculations" begin
    config_file = joinpath(exp_dir, "test", "amip_test.yml")
    cs = CoupledSimulation(config_file)

    # Step twice to get non-zero wind and humidity in the atmosphere so fluxes are non-zero
    step!(cs)
    step!(cs)

    (; atmos_sim, land_sim) = cs.model_sims
    boundary_space = Interfacer.boundary_space(cs)
    land_fraction = Interfacer.get_field(land_sim, Val(:area_fraction))

    # Check that the turbulent fluxes stored in the atmosphere match the land model,
    # comparing only where land fraction is exactly 1.
    # The atmosphere stores fluxes in sfc_conditions as ClimaCore geometry objects:
    #   ρ_flux_h_tot = (F_lh + F_sh) * surface_normal   (C3 vector, 1 component)
    #   ρ_flux_q_tot = F_turb_moisture * surface_normal  (C3 vector, 1 component)
    #   ρ_flux_uₕ   = surface_normal ⊗ C12(ρτxz, ρτyz)  (C3⊗C12 tensor, 2 components)
    atmos_cache = atmos_sim.integrator.p.precomputed.sfc_conditions

    snow_cover = land_sim.integrator.p.snow.snow_cover_fraction
    soil_flux = land_sim.integrator.p.soil.turbulent_fluxes
    snow_flux = land_sim.integrator.p.snow.turbulent_fluxes
    canopy_flux = land_sim.integrator.p.canopy.turbulent_fluxes

    # Atmos geometry objects needed to convert land scalars into atmos flux format
    atmos_surface_space = get_surface_space(atmos_sim)
    Y = atmos_sim.integrator.u
    surface_geometry =
        CC.Fields.level(CC.Fields.local_geometry_field(Y.f), CC.Utilities.half)
    surface_normal = @. CA.C3(CA.unit_basis_vector_data(CA.C3, surface_geometry))
    vec_ct12_ct1 = @. CA.CT12(
        CA.CT1(CA.unit_basis_vector_data(CA.CT1, surface_geometry)),
        surface_geometry,
    )
    vec_ct12_ct2 = @. CA.CT12(
        CA.CT2(CA.unit_basis_vector_data(CA.CT2, surface_geometry)),
        surface_geometry,
    )

    # F_lh + F_sh: compare the single C3 component of (scalar * surface_normal) on both sides
    land_energy_atmos = Interfacer.remap(
        atmos_surface_space,
        canopy_flux.lhf .+ soil_flux.lhf .* (1 .- snow_cover) .+
        snow_cover .* snow_flux.lhf .+ canopy_flux.shf .+
        soil_flux.shf .* (1 .- snow_cover) .+ snow_cover .* snow_flux.shf,
    )
    land_ρ_flux_h_tot = land_energy_atmos .* surface_normal
    err_energy = @. atmos_cache.ρ_flux_h_tot.components.data.:1 -
       land_ρ_flux_h_tot.components.data.:1
    err_energy = @. ifelse(land_fraction < 1, zero(err_energy), err_energy)
    @test maximum(abs.(err_energy)) < 0.1

    # F_turb_moisture (atmosphere stores as scalar * surface_normal)
    ρ_liq = CL.Parameters.ρ_cloud_liq(land_sim.model.soil.parameters.earth_param_set)
    land_moisture_atmos = Interfacer.remap(
        atmos_surface_space,
        (
            canopy_flux.vapor_flux .+
            (soil_flux.vapor_flux_liq .+ soil_flux.vapor_flux_ice) .* (1 .- snow_cover) .+ snow_cover .* snow_flux.vapor_flux
        ) .* ρ_liq,
    )
    land_ρ_flux_q_tot = land_moisture_atmos .* surface_normal
    err_moisture = @. atmos_cache.ρ_flux_q_tot.components.data.:1 -
       land_ρ_flux_q_tot.components.data.:1
    err_moisture = @. ifelse(land_fraction < 1, zero(err_moisture), err_moisture)
    @test maximum(abs.(err_moisture)) < 1e-10

    # F_turb_ρτxz and F_turb_ρτyz: compare the two C3⊗C12 tensor components as scalars
    land_ρτxz_atmos = Interfacer.remap(
        atmos_surface_space,
        soil_flux.ρτxz .* (1 .- snow_cover) .+ snow_cover .* snow_flux.ρτxz,
    )
    land_ρτyz_atmos = Interfacer.remap(
        atmos_surface_space,
        soil_flux.ρτyz .* (1 .- snow_cover) .+ snow_cover .* snow_flux.ρτyz,
    )
    land_ρ_flux_uₕ = (
        surface_normal .⊗
        CA.C12.(
            land_ρτxz_atmos .* vec_ct12_ct1 .+ land_ρτyz_atmos .* vec_ct12_ct2,
            surface_geometry,
        )
    )

    err_ρτ_1 =
        @. atmos_cache.ρ_flux_uₕ.components.data.:1 - land_ρ_flux_uₕ.components.data.:1
    err_ρτ_1 = @. ifelse(land_fraction < 1, zero(err_ρτ_1), err_ρτ_1)
    @test maximum(abs.(err_ρτ_1)) < 1e-10

    err_ρτ_2 =
        @. atmos_cache.ρ_flux_uₕ.components.data.:2 - land_ρ_flux_uₕ.components.data.:2
    err_ρτ_2 = @. ifelse(land_fraction < 1, zero(err_ρτ_2), err_ρτ_2)
    @test maximum(abs.(err_ρτ_2)) < 1e-10
end
