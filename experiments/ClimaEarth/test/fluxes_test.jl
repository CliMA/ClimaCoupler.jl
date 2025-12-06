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
    config_dict = get_coupler_config_dict(config_file)

    # Make sure radiation is computed during the first step
    config_dict["dt_rad"] = config_dict["dt"]

    # Construct coupled simulation and run one coupling step
    cs = CoupledSimulation(config_dict)
    step!(cs)
    boundary_space = Interfacer.boundary_space(cs)
    FT = CC.Spaces.undertype(boundary_space)

    # Unpack component models
    (; atmos_sim, land_sim, ocean_sim, ice_sim) = cs.model_sims

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
    land_flux = Interfacer.remap(land_sim.integrator.p.bucket.R_n, boundary_space)
    @. land_flux = ifelse(land_fraction ≈ 0, zero(FT), land_flux)

    err_land = @. atmos_flux - land_flux
    @. err_land = ifelse(land_fraction ≈ 0, zero(FT), err_land)
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
    @. ice_rad_flux = ifelse(p.area_fraction ≈ 0, zero(FT), ice_rad_flux)
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
    @. ocean_rad_flux = ifelse(ocean_fraction ≈ 0, zero(FT), ocean_rad_flux)

    # Combine component fluxes by area-weighted sum (incl. bucket sign convention):
    combined_fluxes =
        .-land_fraction .* land_flux .+ ice_fraction .* ice_rad_flux .+
        ocean_fraction .* ocean_rad_flux
    err_fluxes = atmos_flux .+ combined_fluxes
    @show "Combined fluxes error: $(maximum(abs.(err_fluxes)))"
    @test maximum(abs.(err_fluxes)) < 8
end
