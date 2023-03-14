struct DiffusiveColumn{FT} <: AbstractComponentModel end

"""
    rhs!(du, u, (parameters, T_sfc), t)

Heat diffusion equation
    dT/dt =  ∇ μ ∇ T
    where
        T  = 280 K              at z = zmax
        dT/dt = - ∇ F_sfc       at z = zmin

We also use this model to calculate and accumulate the downward surface fluxes, F_sfc:
    F_sfc = - λ * (T_sfc - T1)
    d(F_integrated)/dt  = F_sfc
    where
        F_integrated is reset to 0 at the beginning of each coupling cycle
        T1 = atm temperature near the surface (here assumed equal to the first model level)
"""
function get_rhs!(::DiffusiveColumn{FT})
    function rhs!(du, u, cache, t) # see ODE.jl interface

        ## Surface Flux Calculation (coarse bulk formula)
        calculate_flux(T_sfc, T1, p) = -p.λ .* (T_sfc .- T1);

        F_sfc = calculate_flux(cache.T_slab, Fields.level(u.T,1), cache.params)

        ## set BCs
        # bcs_bottom = Operators.SetValue(Geometry.WVector.(F_sfc)) # F_sfc is converted to a Cartesian vector in direction 3 (vertical)
        bcs_bottom = Operators.SetGradient(Geometry.WVector(parent(F_sfc)[1]))
        bcs_top = Operators.SetValue(FT(cache.params.T_top))

        gradc2f = Operators.GradientC2F(top = bcs_top, bottom = bcs_bottom) # Dirichlet BC (center-to-face)
        gradf2c = Operators.DivergenceF2C() # Neumann BC (face-to-center)


        ## tendency calculations
        @. du.T = gradf2c(cache.params.μ * gradc2f(u.T))
        @. du.integrated_flux = - F_sfc # to calculate ∫ F_sfc dt
    end
end

# column init
function init(tag::DiffusiveColumn{FT}; timestepping = (1,1,1)) #; tspan, dt, saveat, space, ocean_mask, stepper = Euler()) where {FT}

    tspan, dt, saveat = timestepping
    # params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))

    # column parameters
    params = (
        zmin = FT(0.0), # height of atm stack bottom [m]
        zmax = FT(1.0), # height of atm stack top [m]
        n = 15,  # number of elements in atm stack
        μ = FT(0.0001), # diffusion coefficient [m^2 / s]
        T_top = FT(280.0), # fixed temperature at the top of the domain [K]
        T_0 = FT(280.0), # initial condition of at temperature (isothermal) [K]
        λ = FT(1e-5), # transfer coefficient
        )

    domain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(params.zmin),
        Geometry.ZPoint{FT}(params.zmax);
        boundary_tags = (:bottom, :top),
    );
    mesh = Meshes.IntervalMesh(domain, nelems = params.n); # struct, allocates face boundaries
    space = Spaces.CenterFiniteDifferenceSpace(mesh); # collection of the above, discretises space into FD and provides coords

    # initialize prognostic variables, either as ClimaCore's Field objects or as Arrays
    T_0 = Fields.ones(FT, space) .* params.T_0; # initiates a spatially uniform atm progostic var

    Y_0 = Fields.FieldVector(; T = T_0, integrated_flux = Fields.level(Fields.zeros(space),1))

    cache = (
        params = params,
        T_slab =Fields.level(Fields.zeros(space),1)
    ) # cache needed for ODEProblem

    problem = OrdinaryDiffEq.ODEProblem(get_rhs!(tag), Y_0, tspan, cache)

    stepper = Euler()
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    ComponentSimulation(tag, Y_init = Y_0, integrator = integrator) # all info needed by the coupler
end







