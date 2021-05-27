using ClimateMachine.Mesh.Grids: _x1, _x2, _x3, VerticalDirection
using ClimateMachine.DGMethods.NumericalFluxes

export CplModel

include("../temp_hooks/overintegration_hook.jl")

struct CplModel{G, D, B, S, TS}
    grid::G
    discretization::D
    boundary::B
    state::S
    odesolver::TS
    nsteps::Int
end

"""
    CplModel(; 
        grid,
        equations,
        nsteps::Int,
        boundary_z = 0.0,
        dt = 1.0,
        timestepper = LSRK54CarpenterKennedy,
        NFfirstorder = RusanovNumericalFlux(),
        NFsecondorder = CentralNumericalFluxSecondOrder(),
        NFgradient = CentralNumericalFluxGradient(),
    ) 

Builds an instance of a coupler test model.  This is a toy model
used for testing and designing coupling machinery. In a full-blown coupled
experiment this model would be replaced by a full compnent model.

-  `grid` the spectral element grid used by this model. 
-  `equations` the Balance Law used by this model.
-  `nsteps` number of component steps to run during each coupling step.
-  `boundary_z` height above or below air-sea interface of the coupled boundary.
-  `dt` component timestep to use on each component step.
-  `timestepper` the ODE solver used to advance the system.
-  `NFfirstorder` numerical flux to use for first order terms.
-  `NFsecondorder` numerical flux to use for second order terms.
-  `NFgradient` numerical flux to use for gradient terms.

Each returned model instance is independent and has its own grid,
balance law, time odesolver and other attributes.  For now the code
keeps some of these things the same for initial testing, including
component timestepper and initial time (both of which need tweaking
to use for real setups).

A real model might have many more flags and/or may wrap the component creation
very differently. Any component should allow itself to set a number of timesteps
to execute with a certain timestep to synchronize with the coupling time scale.

"""
function CplModel(;
    grid,
    equations,
    nsteps::Int,
    boundary_z = nothing,
    dt = 1.0,
    timestepper = LSRK54CarpenterKennedy,
    NFfirstorder = RusanovNumericalFlux(),
    NFsecondorder = CentralNumericalFluxSecondOrder(),
    NFgradient = CentralNumericalFluxGradient(),
    overint_params = nothing,
)

    FT = eltype(grid.vgeo)

    ###
    ### Create a discretization that is the union of the spatial
    ### grid and the equations, plus some numerical flux settings.
    ###
    discretization = DGModel(
        equations,
        grid,
        NFfirstorder,
        NFsecondorder,
        NFgradient,
        direction = EveryDirection(),
    )

    # Specify the coupling boundary
    xc = grid.vgeo[:, _x1:_x1, :]
    yc = grid.vgeo[:, _x2:_x2, :]
    zc = grid.vgeo[:, _x3:_x3, :]
    
    # Default is to set the coupling boundary at z = 0m
    boundary(boundary_z, xc, yc, zc) = isnothing(boundary_z) ? zc .^2 .< eps(FT) : boundary_z( equations.param_set, xc, yc, zc )
    boundary = boundary(boundary_z, xc, yc, zc)
    
    ###
    ### Invoke the spatial ODE initialization functions
    ###
    state = init_ode_state(discretization, FT(0); init_on_cpu = true)

    ###
    ### Additional tendency hooks
    ###
    if overint_params != nothing
        overintegrationorder, polynomialorder = overint_params
        Ns = (polynomialorder.horizontal, polynomialorder.horizontal, polynomialorder.vertical)
        No = (overintegrationorder.horizontal, overintegrationorder.horizontal, overintegrationorder.vertical)
        overintegration_filter!(state, discretization, Ns, No) # only works if Nover > 0
    end 
    function custom_tendency(tendency, x...; kw...)
        discretization(tendency, x...; kw...)
        if overint_params != nothing
            overintegration_filter!(tendency, discretization, Ns, No)
            #"uisng Oi"
        end
    end

    ###
    ### Create a timestepper of the sort needed for this component.
    ### Hard coded here - but can be configurable.
    ###
    odesolver = timestepper(custom_tendency, state, dt = dt, t0 = 0.0)

    ###
    ### Return a CplModel entity that holds all the information
    ### for a component that can be driver from a coupled stepping layer.
    ###
    return CplModel(grid, discretization, boundary, state, odesolver, nsteps)
end

