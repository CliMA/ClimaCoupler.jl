using Test
import StaticArrays
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import ClimaCoupler
import ClimaCoupler.FluxCalculator

# Every trigger package must be imported for ClimaCouplerCMIPExt to load
import Oceananigans,
    ClimaOcean, ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
const CMIPExt = Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt)

const FT = Float64
const sf_params = SurfaceFluxesParameters(FT, UF.BusingerParams)
const thermo_params = SF.Parameters.thermodynamics_params(sf_params)

const σ = FT(5.670374419e-8)
const T_melt = FT(273.15)
const maximum_temperature_step = FT(5)
const height_int = FT(20)
const height_sfc = FT(0)
const z0 = FT(5.8e-5)

"""
    Column(...)

One atmospheric column over sea ice: the atmospheric state at the first interior level,
the downwelling radiation, and the column's conductive resistance `R` and internal
(ice-ocean interface) temperature `T_i`.
"""
Base.@kwdef struct Column
    T_atmos::FT
    u_atmos::FT
    q_atmos::FT = FT(5e-4)
    ρ_atmos::FT = FT(1.3)
    SW_d::FT = FT(0)
    LW_d::FT
    α::FT = FT(0.7)
    ϵ::FT = FT(1)
    R::FT
    T_i::FT = FT(271.35)
end

"""
    step_surface(c::Column, T_guess)

Run one `SurfaceFluxes` solve for the column with the sea-ice skin-temperature callback,
starting from the surface temperature `T_guess`. Returns the coupler's flux NamedTuple,
whose `T_sfc_new` is the updated skin temperature. This mirrors what
`compute_surface_fluxes!(::ClimaSeaIceSimulation, ...)` does for one grid point.
"""
function step_surface(c::Column, T_guess)
    update_T_sfc_cb = CMIPExt.update_T_sfc(
        c.R,
        c.T_i,
        σ,
        c.ϵ,
        c.SW_d,
        c.LW_d,
        c.α,
        T_melt,
        maximum_temperature_step,
    )
    config = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(z0, z0),
        SF.ConstantGustinessSpec(FT(1)),
    )
    uv_int = StaticArrays.SVector(c.u_atmos, zero(FT))
    uv_sfc = StaticArrays.SVector(zero(FT), zero(FT))
    ρ_sfc = SF.surface_density(
        sf_params,
        c.T_atmos,
        c.ρ_atmos,
        T_guess,
        height_int - height_sfc,
        c.q_atmos,
        zero(FT),
        zero(FT),
    )
    q_sfc = TD.q_vap_saturation(thermo_params, T_guess, ρ_sfc, zero(FT), zero(FT))

    return FluxCalculator.get_surface_fluxes(
        sf_params,
        uv_int,
        c.T_atmos,
        c.q_atmos,
        zero(FT),
        zero(FT),
        c.ρ_atmos,
        height_int,
        uv_sfc,
        T_guess,
        q_sfc,
        height_sfc,
        zero(FT),
        config,
        update_T_sfc_cb,
    )
end

"""
    flux_balance_residual(c::Column, fluxes)

Net upward surface flux minus the conductive flux the ice model will draw from the same
skin temperature, in W m⁻². The skin temperature is converged when this vanishes, and the
residual is what ClimaSeaIce turns into spurious surface melt when it does not.
"""
function flux_balance_residual(c::Column, fluxes)
    Tₛ = fluxes.T_sfc_new
    Jᵃ = σ * c.ϵ * Tₛ^4 - (1 - c.α) * c.SW_d - c.ϵ * c.LW_d + fluxes.F_sh + fluxes.F_lh
    return Jᵃ + (Tₛ - c.T_i) / c.R
end

"""
    converge(c::Column, T_guess; n = 40)

Advance the skin temperature the way the coupling loop does: `SurfaceFluxes` applies the
callback once per coupling step, and the result is written back to the ice model and read
again as the next step's guess.
"""
function converge(c::Column, T_guess; n = 40)
    fluxes = step_surface(c, T_guess)
    for _ in 2:n
        fluxes = step_surface(c, fluxes.T_sfc_new)
    end
    return fluxes
end

const columns = (
    "Arctic night, calm" =>
        Column(T_atmos = 250.0, u_atmos = 2.0, LW_d = 180.0, R = 1.0),
    "Arctic night, windy" =>
        Column(T_atmos = 250.0, u_atmos = 15.0, LW_d = 180.0, R = 1.0),
    "thick ice, windy" =>
        Column(T_atmos = 245.0, u_atmos = 15.0, LW_d = 170.0, R = 4.0),
    "thin ice (h = 10 cm)" =>
        Column(T_atmos = 250.0, u_atmos = 10.0, LW_d = 180.0, R = 0.05),
    "new ice (h = 1 cm)" =>
        Column(T_atmos = 250.0, u_atmos = 10.0, LW_d = 180.0, R = 0.005),
    "deep snow" => Column(T_atmos = 240.0, u_atmos = 12.0, LW_d = 160.0, R = 10.0),
    "surface at air temperature" =>
        Column(T_atmos = 260.0, u_atmos = 10.0, LW_d = 220.0, R = 1.0),
)

const melting_column =
    Column(T_atmos = 274.0, u_atmos = 5.0, SW_d = 400.0, LW_d = 300.0, R = 0.5)

@testset "sea-ice skin temperature solver" begin
    @test CMIPExt !== nothing

    @testset "converges to the flux balance: $name" for (name, c) in columns
        fluxes = converge(c, T_melt)
        @test abs(flux_balance_residual(c, fluxes)) < 1  # W m⁻²
        @test 150 < fluxes.T_sfc_new < T_melt
    end

    @testset "independent of the initial guess: $name" for (name, c) in columns
        from_melting = converge(c, T_melt).T_sfc_new
        from_air = converge(c, c.T_atmos).T_sfc_new
        from_cold = converge(c, FT(200)).T_sfc_new
        @test from_melting ≈ from_air atol = 0.05
        @test from_melting ≈ from_cold atol = 0.05
    end

    @testset "no oscillation once converged: $name" for (name, c) in columns
        Tₛ = converge(c, T_melt).T_sfc_new
        for _ in 1:8
            Tₛ⁺ = step_surface(c, Tₛ).T_sfc_new
            @test abs(Tₛ⁺ - Tₛ) < 0.05
            Tₛ = Tₛ⁺
        end
    end

    @testset "step is bounded from any guess" begin
        for (_, c) in columns, T_guess in (FT(150), FT(200), FT(250), T_melt)
            Tₛ⁺ = step_surface(c, T_guess).T_sfc_new
            @test abs(Tₛ⁺ - T_guess) <= maximum_temperature_step + sqrt(eps(FT))
            @test Tₛ⁺ <= T_melt
        end
    end

    @testset "capped at the melting temperature under strong heating" begin
        fluxes = converge(melting_column, T_melt)
        @test fluxes.T_sfc_new ≈ T_melt
        # Heating in excess of what the column can conduct away is the melt energy the
        # ice model absorbs, so the residual must be negative (net downward), not zero.
        @test flux_balance_residual(melting_column, fluxes) < 0
    end
end
