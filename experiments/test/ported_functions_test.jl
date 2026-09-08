using Test
import Oceananigans as OC
import ClimaOcean as CO
import ClimaCoupler

# Every trigger package must be imported for ClimaCouplerCMIPExt to load
import ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
const CMIPExt = Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt)

@testset "ported functions match ClimaOcean" begin
    @test CMIPExt !== nothing
end

# A ported function is a distinct object from its upstream original, so `typeof` never matches even for
# an identical copy. Compare closures structurally instead: same wrapper type, and field-by-field equal
# values, with function-valued fields compared by the values they return rather than by identity.
const closure_sample_points = (
    (0.0, 0.0, 0.0, 0.0),
    (10.0, 30.0, -5.0, 0.0),
    (10.0, -30.0, -50.0, 1e4),
    (-120.0, 45.0, -100.0, 0.0),
    (200.0, -60.0, -500.0, 1e5),
)

const comparison_grid = OC.RectilinearGrid(size = (4, 4, 2), extent = (1, 1, 1))
const discrete_sample_indices = ((1, 1, 1), (2, 3, 1), (4, 4, 2))

function equivalent(ported::Function, upstream::Function)
    nameof(ported) === nameof(upstream) || return false

    if hasmethod(ported, NTuple{4, Float64})
        hasmethod(upstream, NTuple{4, Float64}) || return false
        return all(ported(p...) === upstream(p...) for p in closure_sample_points)
    end

    ℓ = (OC.Center(), OC.Center(), OC.Center())
    discrete_call(f, i, j, k) = f(i, j, k, comparison_grid, ℓ..., nothing, nothing, 10.0)
    if hasmethod(ported, typeof((1, 1, 1, comparison_grid, ℓ..., nothing, nothing, 10.0)))
        hasmethod(
            upstream,
            typeof((1, 1, 1, comparison_grid, ℓ..., nothing, nothing, 10.0)),
        ) || return false
        return all(
            discrete_call(ported, i, j, k) === discrete_call(upstream, i, j, k) for
            (i, j, k) in discrete_sample_indices
        )
    end

    return false
end

function equivalent(ported, upstream)
    typeof(ported).name.wrapper === typeof(upstream).name.wrapper || return false
    names = fieldnames(typeof(ported))
    names === fieldnames(typeof(upstream)) || return false
    isempty(names) && return ported == upstream
    return all(equivalent(getfield(ported, n), getfield(upstream, n)) for n in names)
end

@testset "ocean closures" begin
    for (ported, upstream) in (
        (
            CMIPExt.simplified_ocean_closure(),
            CO.OceanConfigurations.simplified_ocean_closure(),
        ),
        (
            CMIPExt.default_one_degree_closure(),
            CO.OceanConfigurations.default_one_degree_closure(),
        ),
    )
        @test length(ported) == length(upstream)
        @test all(equivalent(p, u) for (p, u) in zip(ported, upstream))
    end

    # The comparison must be able to fail: a closure built with wrong parameters is not equivalent.
    wrong = CO.OceanConfigurations.default_one_degree_closure(
        κ_skew = 999,
        biharmonic_timescale = 1.0,
    )
    reference = CO.OceanConfigurations.default_one_degree_closure()
    @test !all(equivalent(p, u) for (p, u) in zip(wrong, reference))
end

@testset "ocean_simulation" begin
    grid = OC.LatitudeLongitudeGrid(
        OC.CPU();
        size = (180, 90, 8),
        longitude = (0, 360),
        latitude = (-75, 75),
        z = (-5500, 0),
        halo = (5, 5, 4),
    )
    # A seamount rising to -1000 m, so cells are genuinely immersed and the immersed drag is exercised.
    bottom_height = [
        -5500.0 + 4500 * exp(-(((i - 90) / 12)^2 + ((j - 45) / 8)^2)) for
        i in 1:180, j in 1:90
    ]
    grid = OC.ImmersedBoundaryGrid(
        grid,
        OC.GridFittedBottom(bottom_height);
        active_cells_map = true,
    )

    build(f) = f(
        grid;
        free_surface = OC.SplitExplicitFreeSurface(grid; substeps = 70),
        momentum_advection = OC.WENOVectorInvariant(order = 5),
        tracer_advection = OC.WENO(order = 5),
        closure = CMIPExt.simplified_ocean_closure(),
    )

    ported = build(CMIPExt.ocean_simulation)
    upstream = build(CO.ocean_simulation)

    for sim in (ported, upstream)
        OC.set!(
            sim.model,
            u = (λ, φ, z) -> 0.2 * sin(deg2rad(2λ)) * cos(deg2rad(φ)),
            v = (λ, φ, z) -> 0.1 * cos(deg2rad(3λ)),
            T = (λ, φ, z) -> 10 + 1e-3z + 2 * cos(deg2rad(φ)),
            S = 35,
        )
        sim.Δt = 60
        OC.time_step!(sim)
        OC.time_step!(sim)
    end

    u_ported = OC.interior(ported.model.velocities.u)
    @test all(isfinite, u_ported)
    @test maximum(abs, u_ported) > 1e-3

    for name in (:T, :S)
        @test OC.interior(ported.model.tracers[name]) ≈
              OC.interior(upstream.model.tracers[name])
    end
    @test OC.interior(ported.model.velocities.u) ≈ OC.interior(upstream.model.velocities.u)
    @test OC.interior(ported.model.velocities.v) ≈ OC.interior(upstream.model.velocities.v)
    # NumericalEarth wraps the tracer top fluxes in a `MultipleFluxes` discrete boundary function, so
    # `surface_flux` has to unwrap that as well as the plain `FluxBoundaryCondition(field)` form.
    for sim in (ported, upstream)
        for f in (
            sim.model.velocities.u,
            sim.model.velocities.v,
            sim.model.tracers.T,
            sim.model.tracers.S,
        )
            @test CMIPExt.surface_flux(f) isa OC.Field
        end
    end
end

function sea_ice_ocean_test_setup()
    grid = OC.LatitudeLongitudeGrid(
        OC.CPU();
        size = (8, 8, 4),
        longitude = (0, 20),
        latitude = (60, 75),
        z = (-200, 0),
        halo = (7, 7, 7),
    )

    build_ocean() = CO.ocean_simulation(grid; closure = CMIPExt.simplified_ocean_closure())

    ocean_a = build_ocean()
    ocean_b = build_ocean()

    for ocean in (ocean_a, ocean_b)
        OC.set!(
            ocean.model,
            u = (λ, φ, z) -> 0.1 * sin(deg2rad(4λ)),
            v = (λ, φ, z) -> 0.05 * cos(deg2rad(3λ)),
            T = (λ, φ, z) -> -2.5 + 0.02 * (φ - 60) + 5e-3z,
            S = (λ, φ, z) -> 33 + 0.05 * (φ - 60),
        )
    end

    ice = CO.SeaIces.sea_ice_simulation(grid, ocean_a; Δt = 300.0)
    OC.set!(ice.model.ice_concentration, (λ, φ) -> 0.3 + 0.02 * (φ - 60))
    OC.set!(ice.model.ice_thickness, (λ, φ) -> 0.5 + 0.05 * (φ - 60))
    OC.set!(ice.model.velocities.u, (λ, φ) -> 0.05 * sin(deg2rad(2λ)))
    OC.set!(ice.model.velocities.v, (λ, φ) -> 0.03 * cos(deg2rad(2φ)))

    ocean_properties = (; reference_density = 1020, heat_capacity = 3991)

    return grid, ocean_a, ocean_b, ice, ocean_properties
end

# Mirrors the `ocean_ice_interface` NamedTuple assembled in `ClimaSeaIceSimulation`.
function sea_ice_ocean_interface(grid, flux_formulation)
    fluxes = (
        interface_heat = OC.Field{OC.Center, OC.Center, Nothing}(grid),
        frazil_heat = OC.Field{OC.Center, OC.Center, Nothing}(grid),
        salt = OC.Field{OC.Center, OC.Center, Nothing}(grid),
        x_momentum = OC.Field{OC.Face, OC.Center, Nothing}(grid),
        y_momentum = OC.Field{OC.Center, OC.Face, Nothing}(grid),
    )

    return (;
        fluxes,
        flux_formulation,
        temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid),
        salinity = OC.Field{OC.Center, OC.Center, Nothing}(grid),
    )
end

@testset "sea ice/ocean fluxes ($(friction_velocity_name))" for (
    friction_velocity_name,
    ported_friction_velocity,
    upstream_friction_velocity,
) in (
    ("constant friction velocity", 0.002, 0.002),
    (
        "momentum based friction velocity",
        CMIPExt.MomentumBasedFrictionVelocity(),
        CO.InterfaceComputations.MomentumBasedFrictionVelocity(),
    ),
)

    grid, ocean_a, ocean_b, ice, ocean_properties = sea_ice_ocean_test_setup()

    interface_ported = sea_ice_ocean_interface(
        grid,
        CMIPExt.ThreeEquationHeatFlux(ice; friction_velocity = ported_friction_velocity),
    )
    interface_upstream = sea_ice_ocean_interface(
        grid,
        CO.InterfaceComputations.ThreeEquationHeatFlux(
            ice;
            friction_velocity = upstream_friction_velocity,
        ),
    )

    CMIPExt.compute_sea_ice_ocean_fluxes!(
        interface_ported,
        ocean_a,
        ice,
        ocean_properties;
        Δt = ice.Δt,
    )
    CO.InterfaceComputations.compute_sea_ice_ocean_fluxes!(
        interface_upstream,
        ocean_b,
        ice,
        ocean_properties,
    )

    for flux in (:frazil_heat, :interface_heat, :x_momentum)
        @test any(!iszero, OC.interior(interface_upstream.fluxes[flux]))
        @test all(isfinite, OC.interior(interface_upstream.fluxes[flux]))
    end

    for flux in (:frazil_heat, :interface_heat, :salt, :x_momentum, :y_momentum)
        @test OC.interior(interface_ported.fluxes[flux]) ≈
              OC.interior(interface_upstream.fluxes[flux])
    end
    @test OC.interior(interface_ported.temperature) ≈
          OC.interior(interface_upstream.temperature)
    @test OC.interior(interface_ported.salinity) ≈ OC.interior(interface_upstream.salinity)
    @test OC.interior(ocean_a.model.tracers.T) ≈ OC.interior(ocean_b.model.tracers.T)
end

@testset "above_freezing_ocean_temperature!" begin
    grid, ocean_a, ocean_b, ice, _ = sea_ice_ocean_test_setup()

    # The clamp is only meaningful if some of the initial column is below the freezing point.
    liquidus = ice.model.phase_transitions.liquidus
    T_before = copy(OC.interior(ocean_a.model.tracers.T))
    S_before = OC.interior(ocean_a.model.tracers.S)
    Tₘ = ClimaSeaIce.SeaIceThermodynamics.melting_temperature.(Ref(liquidus), S_before)
    @test any(T_before .< Tₘ)

    CMIPExt.above_freezing_ocean_temperature!(ocean_a, grid, ice)
    CO.EarthSystemModels.above_freezing_ocean_temperature!(ocean_b, grid, ice)

    @test OC.interior(ocean_a.model.tracers.T) ≈ OC.interior(ocean_b.model.tracers.T)
    @test !(OC.interior(ocean_a.model.tracers.T) ≈ T_before)
    @test CMIPExt.above_freezing_ocean_temperature!(ocean_a, grid, nothing) === nothing
end

@testset "sea_ice_simulation" begin
    grid, ocean, _, _, _ = sea_ice_ocean_test_setup()

    ported = CMIPExt.sea_ice_simulation(grid, ocean)
    upstream = CO.SeaIces.sea_ice_simulation(grid, ocean)

    for sim in (ported, upstream)
        OC.set!(sim.model.ice_concentration, (λ, φ) -> 0.3 + 0.02 * (φ - 60))
        OC.set!(sim.model.ice_thickness, (λ, φ) -> 0.5 + 0.05 * (φ - 60))
        sim.Δt = 60
        OC.time_step!(sim)
        OC.time_step!(sim)
    end

    h_ported = OC.interior(ported.model.ice_thickness)
    u_ported = OC.interior(ported.model.velocities.u)
    @test any(!iszero, h_ported)
    @test all(isfinite, h_ported)
    @test any(!iszero, u_ported)
    @test maximum(abs, u_ported) > 1e-4

    @test OC.interior(ported.model.ice_thickness) ≈
          OC.interior(upstream.model.ice_thickness)
    @test OC.interior(ported.model.ice_concentration) ≈
          OC.interior(upstream.model.ice_concentration)
    @test OC.interior(ported.model.velocities.u) ≈ OC.interior(upstream.model.velocities.u)
    @test OC.interior(ported.model.velocities.v) ≈ OC.interior(upstream.model.velocities.v)
end
