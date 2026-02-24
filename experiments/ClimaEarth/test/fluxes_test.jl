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

# Use the AMIP setup helpers to construct a coupled simulation
include(joinpath("..", "setup_run.jl"))

@testset "surface radiative flux consistency (AMIP + bucket land)" begin
    # Build AMIP configuration used in CI by default
    config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml")
    config_dict = Input.get_coupler_config_dict(config_file)

    # Make sure radiation is computed during the first step
    config_dict["dt_rad"] = config_dict["dt"]

    config_dict["land_model"] = "integrated"

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
        atmos_sim.integrator.p.radiation.rrtmgp_model.direct_sw_surface_albedo,
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
        atmos_sim.integrator.p.radiation.rrtmgp_model.surface_emissivity,
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
