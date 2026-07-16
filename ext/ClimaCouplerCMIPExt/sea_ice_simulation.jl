# Ported from NumericalEarth v0.5.8, src/SeaIces/sea_ice_simulation.jl

import ClimaSeaIce.SeaIceThermodynamics: ConductiveFlux, IceWaterThermalEquilibrium
import Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper

ocean_reference_density(ocean::OceananigansModelSimulations, FT) = convert(FT, ocean.model.buoyancy.formulation.equation_of_state.reference_density)
ocean_reference_density(::Nothing, FT) = convert(FT, 1026.0)

function default_snow_thermodynamics(grid)
    FT = eltype(grid)
    snow_conductivity = FT(0.31)
    snow_surface_temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    top_heat_boundary_condition = CSI.PrescribedTemperature(snow_surface_temperature.data)
    return CSI.snow_slab_thermodynamics(grid; conductivity = snow_conductivity, top_heat_boundary_condition)
end

"""
    sea_ice_simulation(grid, ocean=nothing;
                       Δt = 5minutes,
                       ice_salinity = 4, # psu
                       advection = nothing,
                       tracers = (),
                       ice_heat_capacity = 2100, # J kg⁻¹ K⁻¹
                       ice_consolidation_thickness = 0.05, # m
                       sea_ice_density = 900, # kg m⁻³
                       snow_density = 330, # kg m⁻³
                       dynamics = sea_ice_dynamics(grid, ocean),
                       bottom_heat_boundary_condition = nothing,
                       top_heat_boundary_condition = nothing,
                       timestepper = :SplitRungeKutta3,
                       phase_transitions = PhaseTransitions(eltype(grid);
                                                            heat_capacity=ice_heat_capacity,
                                                            density=sea_ice_density),
                       conductivity = 2, # W m⁻¹ K⁻¹
                       internal_heat_flux = ConductiveFlux(; conductivity),
                       snow_thermodynamics = default_snow_thermodynamics(grid),
                       clock = nothing)

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
- `Δt`: time step for the sea ice simulation
- `ice_salinity`: salinity of the sea ice (psu)
- `advection`: optional advection scheme for the sea ice model; if `nothing` (default), no advection
               is applied and only thermodynamics evolve the sea ice state
- `tracers`: optional tracers to include in the sea ice model
- `ice_heat_capacity`: heat capacity of the sea ice (J kg⁻¹ K⁻¹)
- `ice_consolidation_thickness`: thickness threshold for sea ice consolidation (m)
- `sea_ice_density`: density of the sea ice (kg m⁻³)
- `snow_density`: density of the snow (kg m⁻³)
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
- `clock`: clock for the sea ice model. Defaults to `nothing`, in which case the model builds its
           own default clock; pass a `Clock` to control the time type (e.g. when coupling)
"""
function sea_ice_simulation(grid, ocean = nothing;
                            Δt = 5minutes,
                            ice_salinity = 4, # psu
                            advection = nothing,
                            tracers = (),
                            ice_heat_capacity = 2100, # J kg⁻¹ K⁻¹
                            ice_consolidation_thickness = 0.05, # m
                            sea_ice_density = 900, # kg m⁻³
                            snow_density = 330, # kg m⁻³
                            dynamics = sea_ice_dynamics(grid, ocean),
                            bottom_heat_boundary_condition = nothing,
                            top_heat_boundary_condition = nothing,
                            timestepper = :SplitRungeKutta3,
                            phase_transitions = CSI.PhaseTransitions(eltype(grid);
                                                                     heat_capacity = ice_heat_capacity,
                                                                     density = sea_ice_density),
                            conductivity = 2, # W m⁻¹ K⁻¹
                            internal_heat_flux = ConductiveFlux(; conductivity),
                            snow_thermodynamics = default_snow_thermodynamics(grid),
                            clock = nothing)

    # Build consistent boundary conditions for the ice model:
    # - bottom -> flux boundary condition
    # - top -> prescribed temperature boundary condition (calculated in the flux computation)

    if isnothing(top_heat_boundary_condition)
        top_surface_temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid)
        top_heat_boundary_condition = CSI.PrescribedTemperature(top_surface_temperature.data)
    end

    if isnothing(bottom_heat_boundary_condition)
        if isnothing(ocean)
            surface_ocean_salinity = 0
        else
            kᴺ = size(ocean.model.grid, 3)
            surface_ocean_salinity = OC.interior(ocean.model.tracers.S, :, :, kᴺ:kᴺ)
        end
        bottom_heat_boundary_condition = IceWaterThermalEquilibrium(surface_ocean_salinity)
    end

    ice_thermodynamics = CSI.sea_ice_slab_thermodynamics(grid;
                                                         internal_heat_flux,
                                                         top_heat_boundary_condition,
                                                         bottom_heat_boundary_condition)

    bottom_heat_flux = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    top_heat_flux    = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    snowfall         = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Only forward `clock` when supplied so the model keeps its own default otherwise.
    clock_kw = isnothing(clock) ? NamedTuple() : (; clock)

    # Build the sea ice model
    sea_ice_model = CSI.SeaIceModel(grid;
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
                                    top_heat_flux,
                                    clock_kw...)

    verbose = false
    sea_ice = OC.Simulation(sea_ice_model; Δt, verbose)

    return sea_ice
end

default_coriolis(ocean::OC.Simulation) = ocean.model.coriolis
default_coriolis(ocean::Nothing) = OC.HydrostaticSphericalCoriolis(; rotation_rate = default_planet_rotation_rate())

default_solver(grid, ocean) = CSI.SplitExplicitSolver(grid; substeps = 120)

# We assume RK3 has a larger timestep
function default_solver(grid, ocean::OC.Simulation)
    substeps = if ocean.model.timestepper isa SplitRungeKuttaTimeStepper
        240
    else
        120
    end
    return CSI.SplitExplicitSolver(grid; substeps)
end

function sea_ice_dynamics(grid, ocean = nothing;
                          sea_ice_ocean_drag_coefficient = 3.24e-3,
                          rheology = CSI.ElastoViscoPlasticRheology(),
                          coriolis = default_coriolis(ocean),
                          free_drift = nothing,
                          solver = default_solver(grid, ocean))

    # `ocean_surface_velocities(ocean)` is inlined here: the singleton vertical dimension is dropped,
    # unlike the surface salinity slice in `sea_ice_simulation`, which keeps it.
    if isnothing(ocean)
        SSU = OC.Fields.ZeroField()
        SSV = OC.Fields.ZeroField()
    else
        kᴺ = size(ocean.model.grid, 3)
        SSU = view(ocean.model.velocities.u, :, :, kᴺ)
        SSV = view(ocean.model.velocities.v, :, :, kᴺ)
    end

    FT = eltype(grid)
    sea_ice_ocean_drag_coefficient = convert(FT, sea_ice_ocean_drag_coefficient)
    ρₑ = ocean_reference_density(ocean, FT)

    τo  = CSI.SemiImplicitStress(uₑ = SSU, vₑ = SSV, Cᴰ = sea_ice_ocean_drag_coefficient, ρₑ = ρₑ)
    τua = OC.Field{OC.Face, OC.Center, Nothing}(grid)
    τva = OC.Field{OC.Center, OC.Face, Nothing}(grid)

    if isnothing(free_drift)
        free_drift = CSI.StressBalanceFreeDrift((u = τua, v = τva), τo)
    end

    return CSI.SeaIceMomentumEquation(grid;
                                      coriolis,
                                      top_momentum_stress = (u = τua, v = τva),
                                      bottom_momentum_stress = τo,
                                      rheology,
                                      free_drift,
                                      solver)
end
