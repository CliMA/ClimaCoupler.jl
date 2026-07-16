# Ported from NumericalEarth v0.5.8, src/Oceans/{ocean_simulation,Oceans,barotropic_potential_forcing}.jl

import Oceananigans.BoundaryConditions: DefaultBoundaryCondition
import Oceananigans.DistributedComputations: all_reduce
import Oceananigans.Architectures: ReactantState, AbstractArchitecture
import Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ, ∂xᶠᶜᶜ, ∂yᶜᶠᶜ
import Oceananigans.TimeSteppers: AdaptiveVerticallyImplicitDiscretization
import Oceananigans.Units: hours, minutes
import Statistics: mean

const TEOS10EquationOfState = OC.BuoyancyFormulations.SeawaterPolynomials.TEOS10.TEOS10EquationOfState

#####
##### Defaults
#####

@inline default_gravitational_acceleration() = OC.defaults.gravitational_acceleration
@inline default_planet_rotation_rate()       = OC.defaults.planet_rotation_rate

struct Default{V}
    value::V
end

"""
    default_or_override(default::Default, alternative_default=default.value) = alternative_default
    default_or_override(override, alternative_default) = override

Either return `default.value`, an `alternative_default`, or an `override`.

The purpose of this function is to help define constructors with "configuration-dependent" defaults. For example, the default bottom drag
should be 0 for a single column model, but 0.003 for a global model. We therefore need a way to specify both the "normal" default 0.003 as
well as the "alternative default" 0, all while respecting user input and changing this to a new value if specified.
"""
default_or_override(default::Default, possibly_alternative_default = default.value) = possibly_alternative_default
default_or_override(override, alternative_default = nothing) = override

#####
##### Barotropic potential forcing
#####

struct XDirection end
struct YDirection end

struct BarotropicPotentialForcing{D, P}
    direction::D
    potential::P
end

Adapt.adapt_structure(to, bpf::BarotropicPotentialForcing) =
    BarotropicPotentialForcing(Adapt.adapt(to, bpf.direction), Adapt.adapt(to, bpf.potential))

const XDirectionBPF = BarotropicPotentialForcing{<:XDirection}
const YDirectionBPF = BarotropicPotentialForcing{<:YDirection}

@inline (bpf::XDirectionBPF)(i, j, k, grid, clock, fields) = -∂xᶠᶜᶜ(i, j, k, grid, bpf.potential)
@inline (bpf::YDirectionBPF)(i, j, k, grid, clock, fields) = -∂yᶜᶠᶜ(i, j, k, grid, bpf.potential)

#####
##### Utilities
#####

keep_user_boundary_condition(user, default) = user isa DefaultBoundaryCondition ? default : user

merge_boundary_conditions(user, default) = user

"""
    merge_boundary_conditions(user::FieldBoundaryConditions, default::FieldBoundaryConditions)

Merge `user` and `default` boundary conditions side-by-side: every side the user left unspecified (a `DefaultBoundaryCondition`) inherits
the corresponding default side. This allows users to prescribe, for example, only the lateral boundary conditions of a field while retaining
the default surface fluxes, bottom drag, and immersed boundary condition.
"""
function merge_boundary_conditions(user::OC.FieldBoundaryConditions, default::OC.FieldBoundaryConditions)
    return OC.FieldBoundaryConditions(keep_user_boundary_condition(user.west, default.west),
                                      keep_user_boundary_condition(user.east, default.east),
                                      keep_user_boundary_condition(user.south, default.south),
                                      keep_user_boundary_condition(user.north, default.north),
                                      keep_user_boundary_condition(user.bottom, default.bottom),
                                      keep_user_boundary_condition(user.top, default.top),
                                      keep_user_boundary_condition(user.immersed, default.immersed))
end

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds -μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds -μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

# Keep a constant linear drag parameter independent on vertical level
@inline u_immersed_bottom_drag(i, j, k, grid, clock, Φ, μ) = @inbounds -μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, Φ, μ) = @inbounds -μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

@inline build_top_bc(flux_field, ::Nothing) = OC.FluxBoundaryCondition(flux_field)

default_free_surface(grid) = OC.SplitExplicitFreeSurface(grid; cfl = 0.7)

estimate_maximum_Δt(grid::OC.RectilinearGrid) = 30minutes # ?

function estimate_maximum_Δt(grid)
    arch = OC.Architectures.architecture(grid)
    Δx = mean(OC.xspacings(grid))
    Δy = mean(OC.yspacings(grid))
    Δθ = rad2deg(mean([Δx, Δy])) / grid.radius

    # The maximum Δt is roughly 1hours * Δθ, giving:
    # - 60 minutes for a 1 degree ocean
    # - 30 minutes for a 0.5 degree ocean
    # - 15 minutes for a 1/4 degree ocean
    # - 7.5 minutes for a 1/8 degree ocean
    # - 3.75 minutes for a 1/16 degree ocean
    # - 1.875 minutes for a 1/32 degree ocean

    # We set the maximum Δt to 1 hour
    Δt = min(1hours, 1hours * Δθ)

    return all_reduce(min, Δt, arch)
end

const TripolarOfSomeKind =
    Union{OC.TripolarGrid, OC.ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:OC.TripolarGrid}}

function default_free_surface(grid::TripolarOfSomeKind; fixed_Δt = estimate_maximum_Δt(grid), cfl = 0.7)
    free_surface = OC.SplitExplicitFreeSurface(grid; cfl, fixed_Δt)
    return free_surface
end

default_stop_time(grid, clock) = default_stop_time(architecture(grid), clock)
default_stop_time(::AbstractArchitecture, clock) = clock.time isa Number ? Inf : Dates.DateTime(9999, 12, 31, 23, 59, 59)
default_stop_time(::ReactantState, clock) = nothing

hasclosure(closure, ClosureType) = closure isa ClosureType
hasclosure(closure_tuple::Tuple, ClosureType) = any(hasclosure(c, ClosureType) for c in closure_tuple)

const OceananigansModelSimulations =
    Union{OC.Simulation{<:OC.HydrostaticFreeSurfaceModel}, OC.Simulation{<:OC.NonhydrostaticModel}}

"""
    ocean_simulation(grid;
                     Δt = estimate_maximum_Δt(grid),
                     closure = default_ocean_closure(),
                     tracers = (:T, :S),
                     free_surface = default_free_surface(grid),
                     reference_density = 1020,
                     rotation_rate = default_planet_rotation_rate(),
                     gravitational_acceleration = default_gravitational_acceleration(),
                     bottom_drag_coefficient = Default(0.003),
                     forcing = NamedTuple(),
                     additional_surface_fluxes = NamedTuple(),
                     biogeochemistry = nothing,
                     timestepper = :SplitRungeKutta3,
                     coriolis = Default(HydrostaticSphericalCoriolis(; rotation_rate)),
                     momentum_advection = WENOVectorInvariant(time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)),
                     tracer_advection = WENO(order=7, time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)),
                     equation_of_state = TEOS10EquationOfState(; reference_density),
                     boundary_conditions::NamedTuple = NamedTuple(),
                     radiative_forcing = nothing,
                     clock = nothing,
                     warn = true,
                     verbose = false)

Construct and return a hydrostatic ocean simulation tailored to `grid`.

This function assembles an Oceananigans's `HydrostaticFreeSurfaceModel` with physically
consistent defaults for advection, closures, the equation of state, surface fluxes, Coriolis,
barotropic pressure–gradient forcing, boundary conditions, and optional biogeochemistry.
It then wraps the model into an Oceananigans's `Simulation` with the specified timestepping options.

## Behaviour and automatic configuration

### Coriolis
- On spherical grids, an `Oceananigans.Coriolis.HydrostaticSphericalCoriolis` object
  is used by default.
- On rectilinear grids, Coriolis force is disabled unless explicitly provided.

### Single-column grids (`grid.Nx == 1 && grid.Ny == 1`)
- Advection is turned off (`momentum_advection = nothing`, `tracer_advection = nothing`).
- Users may override `bottom_drag_coefficient`, but its default is `0`.
- Immersed boundaries are ignored.

### Bottom drag and immersed boundaries
For multi-column grids:
- Quadratic bottom drag is automatically applied to both `u` and `v`.
- Immersed-boundary bottom drag conditions are constructed for both velocity components.
- Barotropic potential forcings for `u` and `v` are also added automatically, and
  user forcing tuples (e.g. `forcing = (u = ..., v = ...)`) are appended if provided.

### Tracers and closures
- `tracers` defaults to `(:T, :S)`.
- If the closure requires turbulent kinetic energy (e.g. `CATKEVerticalDiffusivity`),
  the turbulent kinetic energy `:e` tracer is automatically added while its advection is disabled.

### Boundary conditions
Default boundary conditions are constructed for `u`, `v`, `T`, and `S`, including
surface fluxes and bottom drag. User-provided boundary conditions override the
defaults on a per-field basis.

## Keyword Arguments

- `Δt`: Timestep used by the `Simulation`. Defaults to the maximum stable timestep estimated from the `grid`.
- `clock`: Clock for the underlying model. Defaults to `nothing`, in which case the
  model builds its own default clock. Pass a `Clock` (e.g. `Clock{Float64}(time=0)` or
  a `DateTime`-based clock) to control the time type, for instance when coupling.
- `stop_time`: the end time of the ocean simulation. Defaults to `Inf`.
- `closure`: A turbulence or mixing closure. Defaults to `default_ocean_closure()`.
- `tracers`: Tuple of tracer names. Defaults to `(:T, :S)`.
- `free_surface`: Free–surface solver. Defaults to `default_free_surface(grid)`.
- `reference_density`: Reference seawater density used by the equation of state.
- `rotation_rate`: Planetary rotation rate used for Coriolis forcing.
- `gravitational_acceleration`: Gravitational acceleration, passed to buoyancy.
- `bottom_drag_coefficient`: Bottom drag coefficient. May be a `Default` wrapper.
- `forcing`: Named tuple of additional forcing(s) for individual fields.
- `additional_surface_fluxes`: Named tuple of additional top boundary flux conditions for any field.
- `biogeochemistry`: A biogeochemical model or `nothing`.
- `timestepper`: Time-stepping scheme; options are `:SplitRungeKutta3` (default), or `:QuasiAdamsBashforth2`.
- `coriolis`: Coriolis object or `Default(...)` wrapper.
- `momentum_advection`: Momentum advection scheme.
- `tracer_advection`: Tracer advection scheme or named tuple of schemes.
- `equation_of_state`: Equation of state object. Defaults to TEOS-10 (`TEOS10EquationOfState`).
- `boundary_conditions`: User-supplied boundary conditions; merged with defaults.
- `radiative_forcing`: Additional temperature forcing; merged into `forcing`. Defaults to `nothing`.
- `warn`: If `true`, warnings are emitted for potentially unintended setups.
- `verbose`: If `true`, prints additional setup information.
"""
function ocean_simulation(grid;
                          Δt = estimate_maximum_Δt(grid),
                          clock = Clock(grid),
                          stop_time = default_stop_time(grid, clock),
                          closure = default_ocean_closure(),
                          tracers = (:T, :S),
                          free_surface = default_free_surface(grid),
                          reference_density = 1020,
                          rotation_rate = default_planet_rotation_rate(),
                          gravitational_acceleration = default_gravitational_acceleration(),
                          bottom_drag_coefficient = Default(0.003),
                          forcing = NamedTuple(),
                          additional_surface_fluxes = NamedTuple(),
                          biogeochemistry = nothing,
                          timestepper = :SplitRungeKutta3,
                          coriolis = Default(OC.HydrostaticSphericalCoriolis(; rotation_rate)),
                          momentum_advection = OC.WENOVectorInvariant(time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5)),
                          tracer_advection = OC.WENO(order = 7, time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5)),
                          equation_of_state = TEOS10EquationOfState(; reference_density),
                          boundary_conditions::NamedTuple = NamedTuple(),
                          radiative_forcing = nothing,
                          warn = true,
                          verbose = false)

    FT = eltype(grid)

    if grid isa OC.RectilinearGrid # turn off Coriolis unless user-supplied
        coriolis = default_or_override(coriolis, nothing)
    else
        coriolis = default_or_override(coriolis)
    end

    # Detect whether we are on a single column grid
    Nx, Ny, _ = size(grid)
    single_column_simulation = Nx == 1 && Ny == 1

    if single_column_simulation
        # Let users put a bottom drag if they want
        bottom_drag_coefficient = default_or_override(bottom_drag_coefficient, zero(grid))

        # Don't let users use advection in a single column model
        tracer_advection = nothing
        momentum_advection = nothing

        # No immersed boundaries in a single column grid
        u_immersed_bc = DefaultBoundaryCondition()
        v_immersed_bc = DefaultBoundaryCondition()
    else
        if warn && !(grid isa OC.ImmersedBoundaryGrid) && verbose
            msg = """Are you totally, 100% sure that you want to build a simulation on

                   $(summary(grid))

                   rather than on an ImmersedBoundaryGrid?
                   """
            @warn msg
        end

        bottom_drag_coefficient = default_or_override(bottom_drag_coefficient)

        u_immersed_drag = OC.FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
        v_immersed_drag = OC.FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)

        u_immersed_bc = OC.ImmersedBoundaries.ImmersedBoundaryCondition(bottom = u_immersed_drag)
        v_immersed_bc = OC.ImmersedBoundaries.ImmersedBoundaryCondition(bottom = v_immersed_drag)

        # Forcing for u, v
        barotropic_potential = OC.Field{OC.Center, OC.Center, Nothing}(grid)
        u_forcing = BarotropicPotentialForcing(XDirection(), barotropic_potential)
        v_forcing = BarotropicPotentialForcing(YDirection(), barotropic_potential)

        :u ∈ keys(forcing) && (u_forcing = (u_forcing, forcing[:u]))
        :v ∈ keys(forcing) && (v_forcing = (v_forcing, forcing[:v]))
        forcing = merge(forcing, (u = u_forcing, v = v_forcing))
    end

    if !isnothing(radiative_forcing)
        if :T ∈ keys(forcing)
            T_forcing = (forcing.T, radiative_forcing)
        else
            T_forcing = radiative_forcing
        end
        forcing = merge(forcing, (; T = T_forcing))
    end

    bottom_drag_coefficient = convert(FT, bottom_drag_coefficient)

    # Set up boundary conditions using Field
    top_zonal_momentum_flux = τˣ = OC.Field{OC.Face, OC.Center, Nothing}(grid)
    top_meridional_momentum_flux = τʸ = OC.Field{OC.Center, OC.Face, Nothing}(grid)
    top_ocean_heat_flux = Jᵀ = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    top_salt_flux = Jˢ = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Merge user-supplied additional fluxes with defaults
    default_additional_fluxes = (u = nothing, v = nothing, T = nothing, S = nothing)
    additional = merge(default_additional_fluxes, additional_surface_fluxes)

    # Construct ocean boundary conditions including surface forcing and bottom drag
    u_top_bc = build_top_bc(τˣ, additional.u)
    v_top_bc = build_top_bc(τʸ, additional.v)
    T_top_bc = build_top_bc(Jᵀ, additional.T)
    S_top_bc = build_top_bc(Jˢ, additional.S)

    u_bot_bc = OC.FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    v_bot_bc = OC.FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)

    default_boundary_conditions = (u = OC.FieldBoundaryConditions(top = u_top_bc, bottom = u_bot_bc, immersed = u_immersed_bc),
                                   v = OC.FieldBoundaryConditions(top = v_top_bc, bottom = v_bot_bc, immersed = v_immersed_bc),
                                   T = OC.FieldBoundaryConditions(top = T_top_bc),
                                   S = OC.FieldBoundaryConditions(top = S_top_bc))

    # Merge boundary conditions side-by-side with preference to user
    merged_boundary_conditions = NamedTuple(name => haskey(default_boundary_conditions, name) ?
                                                    merge_boundary_conditions(boundary_conditions[name], default_boundary_conditions[name]) :
                                                    boundary_conditions[name]
                                            for name in keys(boundary_conditions))

    boundary_conditions = merge(default_boundary_conditions, merged_boundary_conditions)
    buoyancy = OC.SeawaterBuoyancy(; gravitational_acceleration, equation_of_state)

    tracer_advection = NamedTuple(name => tracer_advection for name in tracers)

    if hasclosure(closure, CATKEVerticalDiffusivity)
        # Turn off CATKE tracer advection
        tke_advection = (; e = nothing)
        tracer_advection = merge(tracer_advection, tke_advection)
    end

    # Only forward `clock` when supplied so the model keeps its own default otherwise.
    clock_kw = isnothing(clock) ? NamedTuple() : (; clock)

    ocean_model = OC.HydrostaticFreeSurfaceModel(grid;
                                                 buoyancy,
                                                 closure,
                                                 biogeochemistry,
                                                 tracer_advection,
                                                 momentum_advection,
                                                 tracers,
                                                 timestepper,
                                                 free_surface,
                                                 coriolis,
                                                 forcing,
                                                 boundary_conditions,
                                                 clock_kw...)

    ocean = OC.Simulation(ocean_model; Δt, verbose)

    return ocean
end
