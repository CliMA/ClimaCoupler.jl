using Test
import ClimaCoupler
import ClimaCoupler.Interfacer

# Every trigger package must be imported for ClimaCouplerCMIPExt to load
import Adapt,
    ClimaAtmos, ClimaSeaIce, ConservativeRegridding, KernelAbstractions, Oceananigans
import Oceananigans as OC
const CMIPExt = Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt)

const FT = Float64
const σ = 5.670374419e-8
const C_to_K = 273.15

# Open water, two partial covers, and full cover. A `1/ℵ` error is invisible at ℵ = 1.
const concentrations = [0.0 0.25 0.6 1.0]

"""
    bare_ice_simulation()

A `ClimaSeaIceSimulation` around a standalone `sea_ice_simulation`, with `nothing` in the
fields the top heat flux assembly does not read.
"""
function bare_ice_simulation()
    grid = OC.LatitudeLongitudeGrid(
        size = (4, 4),
        longitude = (0, 10),
        latitude = (0, 10),
        topology = (OC.Bounded, OC.Bounded, OC.Flat),
    )
    ice = CMIPExt.sea_ice_simulation(
        grid;
        clock = OC.TimeSteppers.Clock(grid),
        dynamics = nothing,
    )
    return CMIPExt.ClimaSeaIceSimulation(ice, nothing, nothing, nothing, (; σ, C_to_K), 300.0)
end

"""
    setup!(sim, T_sfc_C, radiative_per_ice)

Set the ice concentration and skin temperature, and seed the top heat flux with the
absorbed radiative term as a grid-cell mean. Return the concentration.
"""
function setup!(sim, T_sfc_C, radiative_per_ice)
    ℵ = OC.interior(sim.ice.model.ice_concentration, :, :, 1)
    for j in axes(ℵ, 2)
        ℵ[:, j] .= vec(concentrations)
    end
    OC.interior(CMIPExt.top_thermodynamics(sim).top_surface_temperature, :, :, 1) .= T_sfc_C
    OC.interior(sim.ice.model.external_heat_fluxes.top, :, :, 1) .= ℵ .* radiative_per_ice
    return ℵ
end

@testset "sea-ice top heat flux" begin
    @testset "the per-ice skin flux is recovered by dividing out the concentration" begin
        sim = bare_ice_simulation()
        T_sfc_C, radiative_per_ice = -20.0, -150.0
        F_lh, F_sh = fill(12.0, 4, 4), fill(30.0, 4, 4)
        ℵ = setup!(sim, T_sfc_C, radiative_per_ice)

        CMIPExt.compute_ice_top_heat_flux!(sim, F_lh, F_sh)

        ϵ = FT(Interfacer.get_field(sim, Val(:emissivity)))
        Jᵃ = σ * ϵ * (T_sfc_C + C_to_K)^4 + radiative_per_ice + 12.0 + 30.0
        Qui = OC.interior(sim.ice.model.external_heat_fluxes.top, :, :, 1)
        for i in eachindex(ℵ)
            ℵ[i] > 0 ? (@test Qui[i] / ℵ[i] ≈ Jᵃ) : (@test Qui[i] == 0)
        end
    end
end
