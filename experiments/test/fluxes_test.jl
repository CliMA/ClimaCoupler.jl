#
# Flux consistency test (AMIP + bucket land)
#
# This test sets up an AMIP coupled simulation (atmosphere + bucket land + prescribed
# ocean + prescribed sea ice), advances one coupling step, and verifies that the
# surface radiative flux seen by the atmosphere matches what each surface model
# computes/stores:
# - Atmosphere: uses `sim.integrator.p.radiation.ᶠradiation_flux` (positive downward).
# - Bucket land: compares against `sim.integrator.p.bucket.R_n` on land-dominant cells.
# - Prescribed ocean: skipped — SST is prescribed so radiative fluxes are not computed
#   in the same way.
# - Prescribed sea ice: not stored directly; compute using cache fields and compare
#   to the atmospheric flux on ice-dominant cells.

import Test: @test, @testset
import ClimaCore as CC

include(joinpath("..", "AMIP", "code_loading.jl"))

import ClimaCoupler: FluxCalculator

@testset "surface radiative flux consistency (AMIP + bucket land)" begin
    # Build AMIP configuration used in CI by default
    config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml")
    config_dict = Input.get_coupler_config_dict(config_file)

    # Make sure radiation is computed during the first step
    config_dict["dt_rad"] = config_dict["dt"]

    # Construct coupled simulation and run one coupling step
    cs = CoupledSimulation(config_dict)
    boundary_space = Interfacer.boundary_space(cs)

    # Unpack component models
    (; atmos_sim, land_sim, ocean_sim, ice_sim) = cs.model_sims
    land_fraction = Interfacer.get_field(land_sim, Val(:area_fraction))

    # Check SWD, LWD, albedo, temp, emissivity between atmos and land before flux calculation
    atmos_swd = Interfacer.get_field(boundary_space, atmos_sim, Val(:SW_d))
    land_swd = Interfacer.remap(boundary_space, land_sim.integrator.p.drivers.SW_d)
    err_swd = @. atmos_swd - land_swd
    err_swd = @. ifelse(land_fraction ≈ 0, zero(err_swd), err_swd)
    @test maximum(abs.(err_swd)) < 1e-10

    atmos_lwd = Interfacer.get_field(boundary_space, atmos_sim, Val(:LW_d))
    land_lwd = Interfacer.remap(boundary_space, land_sim.integrator.p.drivers.LW_d)
    err_lwd = @. atmos_lwd - land_lwd
    err_lwd = @. ifelse(land_fraction ≈ 0, zero(err_lwd), err_lwd)
    @test maximum(abs.(err_lwd)) < 1e-10

    atmos_albedo = CC.Fields.array2field(
        CA.RRTMGP.direct_sw_surface_albedo(atmos_sim.integrator.p.radiation.rrtmgp_model),
        boundary_space,
    )
    land_albedo =
        Interfacer.get_field(boundary_space, land_sim, Val(:surface_direct_albedo))
    err_albedo = @. atmos_albedo - land_albedo
    err_albedo = @. ifelse(land_fraction ≈ 0, zero(err_albedo), err_albedo)
    @test maximum(abs.(err_albedo)) < 1e-10

    atmos_temp = Interfacer.remap(
        boundary_space,
        atmos_sim.integrator.p.precomputed.sfc_conditions.T_sfc,
    )
    land_temp = Interfacer.get_field(boundary_space, land_sim, Val(:surface_temperature))
    err_temp = @. atmos_temp - land_temp
    err_temp = @. ifelse(land_fraction ≈ 0, zero(err_temp), err_temp)
    @test maximum(abs.(err_temp)) < 1e-6

    atmos_emissivity = CC.Fields.array2field(
        CA.RRTMGP.surface_emissivity(atmos_sim.integrator.p.radiation.rrtmgp_model),
        boundary_space,
    )
    land_emissivity = Interfacer.get_field(boundary_space, land_sim, Val(:emissivity))
    err_emissivity = @. atmos_emissivity - land_emissivity
    err_emissivity = @. ifelse(land_fraction ≈ 0, zero(err_emissivity), err_emissivity)
    @test maximum(abs.(err_emissivity)) < 1e-10

    step!(cs)
    boundary_space = Interfacer.boundary_space(cs)

    # Atmosphere: radiative flux on the surface interface
    # Convention: positive downward to the surface
    atmos_flux = CC.Spaces.level(
        atmos_sim.integrator.p.radiation.ᶠradiation_flux.components.data.:1,
        CC.Utilities.half,
    )

    # Bucket land: compare to net radiation stored in the bucket cache
    # Convention note: bucket R_n is stored with the same sign (positive upward,
    # see climaland_bucket.jl), so we compare atmos_flux ≈ R_n at land-dominant points.
    p = land_sim.integrator.p
    land_fraction = Interfacer.get_field(land_sim, Val(:area_fraction))
    land_flux = Interfacer.remap(boundary_space, land_sim.integrator.p.bucket.R_n)
    @. land_flux = ifelse(land_fraction ≈ 0, zero(land_flux), land_flux)

    err_land = @. atmos_flux - land_flux
    @. err_land = ifelse(land_fraction ≈ 0, zero(err_land), err_land)
    @show "Bucket flux error: $(maximum(abs.(err_land)))"
    @test maximum(abs.(err_land)) < 5

    # Prescribed ice: radiative fluxes aren't stored; compute from cache and compare
    p = ice_sim.integrator.p
    Y = ice_sim.integrator.u
    FT = eltype(Y)

    # Radiative flux toward surface (positive downward)
    # TODO: get sigma from parameters
    σ = FT(5.67e-8)
    (; k_ice, h, T_base, ρ, c, α, ϵ) = p.params
    ice_rad_flux =
        (1 .- α) .* p.SW_d .+
        ϵ .* (p.LW_d .- σ .* Interfacer.get_field(ice_sim, Val(:surface_temperature)) .^ 4)
    @. ice_rad_flux = ifelse(p.area_fraction ≈ 0, zero(ice_rad_flux), ice_rad_flux)
    ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))

    # Prescribed ocean: SST is prescribed, but for this test we can still compute
    # the radiative flux seen by the ocean surface using the same formula.
    α = Interfacer.get_field(ocean_sim, Val(:surface_direct_albedo))
    ϵ = Interfacer.get_field(ocean_sim, Val(:emissivity))
    ocean_rad_flux =
        (1 .- α) .* cs.fields.SW_d .+
        ϵ .* (
            cs.fields.LW_d .-
            σ .* Interfacer.get_field(ocean_sim, Val(:surface_temperature)) .^ 4
        )
    ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction))
    @. ocean_rad_flux = ifelse(ocean_fraction ≈ 0, zero(ocean_rad_flux), ocean_rad_flux)

    # Combine component fluxes by area-weighted sum (incl. bucket sign convention):
    combined_fluxes =
        .-land_fraction .* land_flux .+ ice_fraction .* ice_rad_flux .+
        ocean_fraction .* ocean_rad_flux

    @info "Combined fluxes: $(combined_fluxes)"
    @info "Atmos flux: $(atmos_flux)"
    err_fluxes = atmos_flux .+ combined_fluxes
    @show "Combined fluxes error: $(maximum(abs.(err_fluxes)))"
    @test maximum(abs.(err_fluxes)) < 10
end

# End-to-end test for the FluxAccumulator setup, exercising both the slow slab
# ocean and the slow bucket land in one CoupledSimulation. `slabplanet` mode
# allocates both a SlabOceanSimulation and a BucketSimulation, so a single
# (expensive-to-construct) `CoupledSimulation` covers both accumulator paths.
# We override timesteps so both surfaces are slow (dt_ocean = dt_land = 3 *
# dt_cpl) and verify that:
#   1) FluxAccumulators are allocated for both surfaces (and only those),
#   2) each surface's turbulent-flux BC is set by the initial
#      `push_ready_accumulators!` call (so neither surface sees zero fluxes
#      at t=0),
#   3) the accumulators collect contributions across coupling steps and are
#      averaged and pushed only when the slow surface is about to step.
@testset "flux accumulation: slow slab ocean and slow bucket land" begin
    config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/slabplanet_default.yml")
    config_dict = Input.get_coupler_config_dict(config_file)

    # Override timesteps so the ocean and land are "slow"
    config_dict["dt_cpl"] = "200secs"
    config_dict["dt_atmos"] = config_dict["dt_seaice"] = config_dict["dt_cpl"]
    config_dict["dt_ocean"] = config_dict["dt_land"] = "600secs"
    haskey(config_dict, "dt") && delete!(config_dict, "dt")
    # Call `parse_component_dts!` to populate `config_dict["component_dt_dict"]`
    Input.parse_component_dts!(config_dict)
    # Two slow-surface steps within the test window:
    config_dict["t_end"] = "1200secs"

    cs = CoupledSimulation(config_dict)

    # 1) Accumulators allocated for both slow surfaces, and only those.
    @test Set(keys(cs.flux_accumulators)) == Set([:ocean_sim, :land_sim])
    @test Interfacer.sim_dt(cs.model_sims.ocean_sim) > Float64(cs.Δt_cpl)
    @test Interfacer.sim_dt(cs.model_sims.land_sim) > Float64(cs.Δt_cpl)

    # 2) After initialization, each slow surface's turbulent-flux BC should be
    # nonzero: the constructor calls `push_ready_accumulators!(...; force = true)`
    # to populate all slow-surface BCs before the run loop.
    F_turb_ocean_after_init =
        copy(parent(cs.model_sims.ocean_sim.integrator.p.F_turb_energy))
    F_lh_bucket_after_init =
        copy(parent(cs.model_sims.land_sim.integrator.p.bucket.turbulent_fluxes.lhf))
    @test any(F_turb_ocean_after_init .!= 0)
    @test any(F_lh_bucket_after_init .!= 0)
    # And both accumulators were reset after the init push.
    @test cs.flux_accumulators.ocean_sim.n_steps[] == 0
    @test cs.flux_accumulators.land_sim.n_steps[] == 0

    # 3) Take coupling steps and check accumulator / surface-BC behavior.
    # After step 1 (t = 200): turbulent_fluxes! adds one contribution to each
    # accumulator. The next coupling time is t=400, and will_step(_, 400) is
    # false for both (slow dt=600 > 400-0), so no push yet.
    step!(cs)
    @test cs.flux_accumulators.ocean_sim.n_steps[] == 1
    @test cs.flux_accumulators.land_sim.n_steps[] == 1
    @test parent(cs.model_sims.ocean_sim.integrator.p.F_turb_energy) ≈
          F_turb_ocean_after_init
    @test parent(cs.model_sims.land_sim.integrator.p.bucket.turbulent_fluxes.lhf) ≈
          F_lh_bucket_after_init

    # After step 2 (t = 400): turbulent_fluxes! adds a second contribution to
    # each accumulator, then checks will_step(_, 600). Since 600 - 0 >= 600
    # for both surfaces, each accumulator averages the two contributions,
    # pushes them to its surface's BC, and resets.
    step!(cs)
    @test cs.flux_accumulators.ocean_sim.n_steps[] == 0
    @test cs.flux_accumulators.land_sim.n_steps[] == 0

    # After step 3 (t = 600): step_model_sims! advances the slow surfaces (their
    # BCs already hold the 2-step average from the end of step 2). Then
    # turbulent_fluxes! adds one new contribution to each accumulator and
    # checks will_step(_, 800): 800 - 600 < 600, so no push yet.
    step!(cs)
    @test cs.flux_accumulators.ocean_sim.n_steps[] == 1
    @test cs.flux_accumulators.land_sim.n_steps[] == 1
end
