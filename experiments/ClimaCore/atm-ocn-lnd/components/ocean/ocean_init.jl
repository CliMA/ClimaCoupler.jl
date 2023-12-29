
using Random
using Printf
using Oceananigans
using Oceananigans.Units: minute, minutes, hour


#####
##### Parameters
#####

"""
    ocean_init(; kwargs...)

Return an `Oceananigans.Simulation` of a column model initialized with
mutable array surface boundary conditions and a linear density stratification.

Arguments
=========

    * Nz: Number of vertical grid points
    * Lz: [m] Vertical extent of domain,
    * f:  [s-1] Coriolis parameter,
    * g:  [m s-2] Gravitational acceleration,
    * T₀: [ᵒC], sea surface temperature,
    * S₀: [psu], sea surface salinity,
    * α: Thermal expansion coefficient
    * β: Haline contraction coefficient
"""
function ocean_init(
    T_sfc_init; # K
    Nz = 64,  # Number of vertical grid points
    Lz = 512, # Vertical extent of domain
    f = 1e-4, # Coriolis parameter
    g = 9.81, # Gravitational acceleration
    α = 2e-4, # Thermal expansion coefficient
    β = 8e-5, # Haline contraction coefficient
)

    # L = Lz
    # topology = (Periodic, Flat, Bounded)
    # grid = RectilinearGrid(topology=topology, size=(1, 1, Nz), halo=(5, 5, 5), extent=(L, L, L))

    #     bottom(x, y) = x > 80kilometers && x < 90kilometers ? 0.0 : -500meters

    # # grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))

    # coriolis = FPlane(f = 1e-4)

    # model = HydrostaticFreeSurfaceModel(grid = grid,
    #                                     coriolis = coriolis,
    #                                     free_surface = free_surface)


    # u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition([0.0]))
    # v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition([0.0]))
    # T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition([0.0]))
    # S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition([0.0]))

    # eos = LinearEquationOfState(α = α, β = β)

    # model = HydrostaticFreeSurfaceModel(
    #     grid = grid,
    #     tracers = (:T, :S, :e),
    #     free_surface = ImplicitFreeSurface(gravitational_acceleration = g),
    #     buoyancy = SeawaterBuoyancy(gravitational_acceleration = g, equation_of_state = eos),
    #     coriolis = FPlane(f = f),
    #     boundary_conditions = (T = T_bcs, S = S_bcs, u = u_bcs, v = v_bcs),
    #     closure = TKEBasedVerticalDiffusivity(),
    # )

    # simulation = Oceananigans.Simulation(model, Δt = 0.02, stop_iteration = 1)

    # # Initialize the ocean state with a linear temperature and salinity stratification
    # α = simulation.model.buoyancy.model.equation_of_state.α
    # β = simulation.model.buoyancy.model.equation_of_state.β
    # Tᵢ(x, y, z) = 16 + α * g * 5e-5 * z
    # Sᵢ(x, y, z) = 35 - β * g * 5e-5 * z
    # set!(simulation.model, T = Tᵢ, S = Sᵢ)

    # # collect data (needs optimising)
    # ocean_data = []
    # data = (T = deepcopy(simulation.model.tracers.T), time = simulation.model.clock.time)
    # push!(ocean_data, data)


    # USE EXISTING EXAMPLE: ocean_wind_mixing_and_convection.jl

    # ## The grid
    #
    # We use 1×1×24 grid points with 2 m grid with varying spacing in the vertical
    # with higher resolution closer to the
    # surface. Here we use a stretching function for the vertical nodes that
    # maintains relatively constant vertical spacing in the mixed layer, which
    # is desirable from a numerical standpoint:

    Nx = Ny = 1     # number of points in each of horizontal directions
    Nz = 24          # number of points in the vertical direction

    Lx = Ly = 64     # (m) domain horizontal extents
    Lz = 32          # (m) domain depth

    refinement = 1.2 # controls spacing near surface (higher means finer spaced)
    stretching = 12  # controls rate of stretching at bottom

    ## Normalized height ranging from 0 to 1
    h(k) = (k - 1) / Nz

    ## Linear near-surface generator
    ζ₀(k) = 1 + (h(k) - 1) / refinement

    ## Bottom-intensified stretching function
    Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

    ## Generating function
    z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

    grid = RectilinearGrid(size = (Nx, Nx, Nz),
                            x = (0, Lx),
                            y = (0, Ly),
                            z = z_faces)

    # ## Buoyancy that depends on temperature and salinity
    #
    # We use the `SeawaterBuoyancy` model with a linear equation of state,

    buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = 2e-4,
                                                                        haline_contraction = 8e-4))

    # ## Boundary conditions
    #
    # We calculate the surface temperature flux associated with surface cooling of
    # 200 W m⁻², reference density `ρₒ`, and heat capacity `cᴾ`,

    Qʰ = 200.0  # W m⁻², surface _heat_ flux
    ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
    cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater

    Qᵀ = Qʰ / (ρₒ * cᴾ) # K m s⁻¹, surface _temperature_ flux

    # Finally, we impose a temperature gradient `dTdz` both initially and at the
    # bottom of the domain, culminating in the boundary conditions on temperature,

    dTdz = 0.01 # K m⁻¹

    T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ),
                                    bottom = GradientBoundaryCondition(dTdz))

    # Note that a positive temperature flux at the surface of the ocean
    # implies cooling. This is because a positive temperature flux implies
    # that temperature is fluxed upwards, out of the ocean.
    #
    # For the velocity field, we imagine a wind blowing over the ocean surface
    # with an average velocity at 10 meters `u₁₀`, and use a drag coefficient `cᴰ`
    # to estimate the kinematic stress (that is, stress divided by density) exerted
    # by the wind on the ocean:

    u₁₀ = 10    # m s⁻¹, average wind velocity 10 meters above the ocean
    cᴰ = 2.5e-3 # dimensionless drag coefficient
    ρₐ = 1.225  # kg m⁻³, average density of air at sea-level

    Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²

    # The boundary conditions on `u` are thus

    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

    # For salinity, `S`, we impose an evaporative flux of the form

    @inline Qˢ(x, y, t, S, evaporation_rate) = - evaporation_rate * S # [salinity unit] m s⁻¹

    # where `S` is salinity. We use an evporation rate of 1 millimeter per hour,

    evaporation_rate = 1e-3 / hour # m s⁻¹

    # We build the `Flux` evaporation `BoundaryCondition` with the function `Qˢ`,
    # indicating that `Qˢ` depends on salinity `S` and passing
    # the parameter `evaporation_rate`,

    evaporation_bc = FluxBoundaryCondition(Qˢ, field_dependencies=:S, parameters=evaporation_rate)

    # The full salinity boundary conditions are

    S_bcs = FieldBoundaryConditions(top=evaporation_bc)

    # ## Model instantiation
    #
    # We fill in the final details of the model here: upwind-biased 5th-order
    # advection for momentum and tracers, 3rd-order Runge-Kutta time-stepping,
    # Coriolis forces, and the `AnisotropicMinimumDissipation` closure
    # for large eddy simulation to model the effect of turbulent motions at
    # scales smaller than the grid scale that we cannot explicitly resolve.

    model = NonhydrostaticModel(; grid, buoyancy,
                                advection = UpwindBiasedFifthOrder(),
                                timestepper = :RungeKutta3,
                                tracers = (:T, :S),
                                coriolis = FPlane(f=1e-4),
                                closure = AnisotropicMinimumDissipation(),
                                boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs))

    # Notes:
    #
    # * To use the Smagorinsky-Lilly turbulence closure (with a constant model coefficient) rather than
    #   `AnisotropicMinimumDissipation`, use `closure = SmagorinskyLilly()` in the model constructor.
    #
    # * To change the architecture to `GPU`, replace `CPU()` with `GPU()` inside the
    #   `grid` constructor.

    # ## Initial conditions
    #
    # Our initial condition for temperature consists of a linear stratification superposed with
    # random noise damped at the walls, while our initial condition for velocity consists
    # only of random noise.

    ## Random noise damped at top and bottom
    Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

    ## Temperature initial condition: a stable density gradient with random noise superposed.
    T0 = T_sfc_init - 273  # C
    Tᵢ(x, y, z) = T0 + dTdz * z + dTdz * model.grid.Lz * 1e-6 * Ξ(z)

    ## Velocity initial condition: random noise scaled by the friction velocity.
    uᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-3 * Ξ(z)

    ## `set!` the `model` fields using functions or constants:
    set!(model, u=uᵢ, w=uᵢ, T=Tᵢ, S=35)

    # ## Setting up a simulation
    #
    # We set-up a simulation with an initial time-step of 10 seconds
    # that stops at 40 minutes, with adaptive time-stepping and progress printing.

    simulation = Simulation(model, Δt=10.0, stop_time=40minutes)


    return simulation
end
