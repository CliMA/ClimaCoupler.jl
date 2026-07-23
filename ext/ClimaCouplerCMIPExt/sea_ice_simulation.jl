import ClimaSeaIce: default_sea_ice_boundary_conditions
import ClimaSeaIce.SeaIceThermodynamics: ConductiveFlux, IceWaterThermalEquilibrium
import ClimaSeaIce.SeaIceDynamics: maybe_extended_grid
import Oceananigans.BoundaryConditions: Zipper

ocean_reference_density(ocean::OceananigansModelSimulations, FT) =
    convert(FT, ocean.model.buoyancy.formulation.equation_of_state.reference_density)
ocean_reference_density(::Nothing, FT) = convert(FT, 1026.0)

@inline u_immersed_side_drag(i, j, k, grid, clock, Φ, μ) =
    @inbounds -μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_side_drag(i, j, k, grid, clock, Φ, μ) =
    @inbounds -μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

function default_snow_thermodynamics(grid)
    FT = eltype(grid)
    snow_conductivity = FT(0.31)
    snow_surface_temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    top_heat_boundary_condition = CSI.PrescribedTemperature(snow_surface_temperature.data)
    return CSI.snow_slab_thermodynamics(
        grid;
        conductivity = snow_conductivity,
        top_heat_boundary_condition,
    )
end

correct_tripolar_bcs(grid, bcs) = bcs

function correct_tripolar_bcs(grid::TripolarGridOfSomeKind, bcs)
    if bcs.north isa OC.BoundaryCondition && bcs.north.classification isa Zipper
        north = OC.BoundaryCondition(bcs.north.classification, -bcs.north.condition)
        bcs = OC.FieldBoundaryConditions(
            bcs.west,
            bcs.east,
            bcs.south,
            north,
            bcs.bottom,
            bcs.top,
            bcs.immersed,
        )
    end
    return bcs
end

"""
    sea_ice_simulation(grid, ocean=nothing;
                       clock = Clock(grid),
                       stop_time = default_stop_time(grid, clock),
                       Δt = 5minutes,
                       ice_salinity = 4, # psu
                       advection = nothing,
                       tracers = (),
                       ice_heat_capacity = 2100, # J kg⁻¹ K⁻¹
                       ice_consolidation_thickness = 0.05, # m
                       sea_ice_density = 900, # kg m⁻³
                       snow_density = 330, # kg m⁻³
                       side_drag_coefficient = 3e-3,
                       dynamics = sea_ice_dynamics(grid, ocean),
                       bottom_heat_boundary_condition = nothing,
                       top_heat_boundary_condition = nothing,
                       timestepper = :SplitRungeKutta3,
                       phase_transitions = PhaseTransitions(eltype(grid);
                                                            heat_capacity=ice_heat_capacity,
                                                            density=sea_ice_density),
                       conductivity = 2, # W m⁻¹ K⁻¹
                       internal_heat_flux = ConductiveFlux(; conductivity),
                       snow_thermodynamics = default_snow_thermodynamics(grid))

Construct a sea ice simulation with the given grid and optional ocean simulation.
The sea ice model is configured with a slab thermodynamics, Elasto-Visco-Plastic rheology,
and a SplitExplicit Runge-Kutta 3rd order time stepper by default. The thermodynamics
include conductive internal heat flux, and the option to specify top and bottom heat
boundary conditions. The dynamics include a semi-implicit ocean stress formulation,
with the option to specify a free drift velocity.

Arguments
=========
- `grid`: the grid on which to build the sea ice model
- `ocean`: optional ocean simulation to provide surface velocities and salinity for the sea ice

Keyword Arguments
=================
- `clock`: Clock for the underlying model. Defaults to `Clock(grid)`, a numeric clock starting at `time = 0`.
  Pass a `DateTime`-based clock to step the simulation in calendar time (e.g. when coupling).
- `stop_time`: Stop time for the simulation. Defaults to `Inf` for numeric clocks, or
  `DateTime(9999, 12, 31, 23, 59, 59)` for `DateTime` clocks. On Reactant architectures it defaults to `nothing`,
  since Reactant does not support `stop_time`.
- `Δt`: time step for the sea ice simulation
- `ice_salinity`: salinity of the sea ice (psu)
- `advection`: optional advection scheme for the sea ice model; if `nothing` (default), no advection
               is applied and only thermodynamics evolve the sea ice state
- `tracers`: optional tracers to include in the sea ice model
- `ice_heat_capacity`: heat capacity of the sea ice (J kg⁻¹ K⁻¹)
- `ice_consolidation_thickness`: thickness threshold for sea ice consolidation (m)
- `sea_ice_density`: density of the sea ice (kg m⁻³)
- `snow_density`: density of the snow (kg m⁻³)
- `side_drag_coefficient`: lateral (immersed-boundary) drag coefficient for the sea ice velocities
- `dynamics`: sea ice dynamics model to use (default is `sea_ice_dynamics(grid, ocean)`)
- `bottom_heat_boundary_condition`: heat boundary condition at the ice-ocean interface (default
                                    is `IceWaterThermalEquilibrium` with ocean surface salinity)
- `top_heat_boundary_condition`: heat boundary condition at the ice-atmosphere interface (default
                                 is a prescribed temperature calculated in the flux computation)
- `timestepper`: time stepper to use for the sea ice model (default is `:SplitRungeKutta3`)
- `phase_transitions`: phase transition properties for the sea ice (default is a `PhaseTransitions`
                       with specified heat capacity and density)
- `conductivity`: thermal conductivity for the internal heat flux (W m⁻¹ K⁻¹)
- `internal_heat_flux`: internal heat flux formulation for the sea ice (default is a
                        `ConductiveFlux` with specified conductivity)
- `snow_thermodynamics`: thermodynamics for the snow layer (default is a slab thermodynamics with
                         specified conductivity and prescribed temperature)
"""
function sea_ice_simulation(
    grid,
    ocean = nothing;
    clock = Clock(grid),
    stop_time = default_stop_time(grid, clock),
    Δt = 5minutes,
    ice_salinity = 4, # psu
    advection = nothing,
    tracers = (),
    ice_heat_capacity = 2100, # J kg⁻¹ K⁻¹
    ice_consolidation_thickness = 0.05, # m
    sea_ice_density = 900, # kg m⁻³
    snow_density = 330, # kg m⁻³
    side_drag_coefficient = 3e-3,
    dynamics = sea_ice_dynamics(grid, ocean),
    bottom_heat_boundary_condition = nothing,
    top_heat_boundary_condition = nothing,
    timestepper = :SplitRungeKutta3,
    phase_transitions = CSI.PhaseTransitions(
        eltype(grid);
        heat_capacity = ice_heat_capacity,
        density = sea_ice_density,
    ),
    conductivity = 2, # W m⁻¹ K⁻¹
    internal_heat_flux = ConductiveFlux(; conductivity),
    snow_thermodynamics = default_snow_thermodynamics(grid),
)

    # Build consistent boundary conditions for the ice model:
    # - bottom -> flux boundary condition
    # - top -> prescribed temperature boundary condition (calculated in the flux computation)

    if isnothing(top_heat_boundary_condition)
        top_surface_temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid)
        top_heat_boundary_condition =
            CSI.PrescribedTemperature(top_surface_temperature.data)
    end

    if isnothing(bottom_heat_boundary_condition)
        surface_ocean_salinity = ocean_surface_salinity(ocean)
        bottom_heat_boundary_condition = IceWaterThermalEquilibrium(surface_ocean_salinity)
    end

    ice_thermodynamics = CSI.sea_ice_slab_thermodynamics(
        grid;
        internal_heat_flux,
        top_heat_boundary_condition,
        bottom_heat_boundary_condition,
    )

    bottom_heat_flux = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    top_heat_flux = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    snowfall = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    side_drag_coefficient = convert(eltype(grid), side_drag_coefficient)
    u_bc = OC.FluxBoundaryCondition(
        u_immersed_side_drag,
        discrete_form = true,
        parameters = side_drag_coefficient,
    )
    v_bc = OC.FluxBoundaryCondition(
        v_immersed_side_drag,
        discrete_form = true,
        parameters = side_drag_coefficient,
    )

    immersed_u_bc =
        OC.ImmersedBoundaries.ImmersedBoundaryCondition(south = u_bc, north = u_bc)
    immersed_v_bc =
        OC.ImmersedBoundaries.ImmersedBoundaryCondition(west = v_bc, east = v_bc)

    u_bcs = correct_tripolar_bcs(
        grid,
        OC.FieldBoundaryConditions(
            grid,
            (OC.Face(), OC.Center(), nothing);
            immersed = immersed_u_bc,
        ),
    )
    v_bcs = correct_tripolar_bcs(
        grid,
        OC.FieldBoundaryConditions(
            grid,
            (OC.Center(), OC.Face(), nothing);
            immersed = immersed_v_bc,
        ),
    )

    # Build the sea ice model
    sea_ice_model = CSI.SeaIceModel(
        grid;
        clock,
        ice_salinity,
        advection,
        tracers,
        ice_consolidation_thickness,
        sea_ice_density,
        snow_density,
        phase_transitions,
        ice_thermodynamics,
        snow_thermodynamics,
        snowfall,
        dynamics,
        timestepper,
        bottom_heat_flux,
        boundary_conditions = (u = u_bcs, v = v_bcs),
        top_heat_flux,
    )

    verbose = false
    sea_ice = OC.Simulation(sea_ice_model; Δt, stop_time, verbose)

    return sea_ice
end

default_coriolis(ocean::OC.Simulation) = ocean.model.coriolis
default_coriolis(ocean::Nothing) =
    OC.HydrostaticSphericalCoriolis(; rotation_rate = default_planet_rotation_rate())

function sea_ice_dynamics(
    grid,
    ocean = nothing;
    sea_ice_ocean_drag_coefficient = 3.24e-3,
    rheology = CSI.ElastoViscoPlasticRheology(),
    coriolis = default_coriolis(ocean),
    free_drift = nothing,
    solver = CSI.SplitExplicitSolver(grid; substeps = 100),
)

    SSU, SSV = ocean_surface_velocities(ocean)

    FT = eltype(grid)
    sea_ice_ocean_drag_coefficient = convert(FT, sea_ice_ocean_drag_coefficient)
    ρₑ = ocean_reference_density(ocean, FT)

    τo = CSI.SemiImplicitStress(
        uₑ = SSU,
        vₑ = SSV,
        Cᴰ = sea_ice_ocean_drag_coefficient,
        ρₑ = ρₑ,
    )

    velocity_grid = maybe_extended_grid(solver, grid)

    τua = OC.Field{OC.Face, OC.Center, Nothing}(
        velocity_grid;
        boundary_conditions = default_sea_ice_boundary_conditions(velocity_grid, :u),
    )
    τva = OC.Field{OC.Center, OC.Face, Nothing}(
        velocity_grid;
        boundary_conditions = default_sea_ice_boundary_conditions(velocity_grid, :v),
    )

    if isnothing(free_drift)
        free_drift = CSI.StressBalanceFreeDrift((u = τua, v = τva), τo)
    end

    return CSI.SeaIceMomentumEquation(
        velocity_grid;
        coriolis,
        top_momentum_stress = (u = τua, v = τva),
        bottom_momentum_stress = τo,
        rheology,
        free_drift,
        solver,
    )
end
