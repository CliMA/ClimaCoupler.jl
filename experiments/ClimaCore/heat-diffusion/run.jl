# # Heat Equation + Slab Tutorial

# In this tutorial, we demonstrate simple sequential coupling
# of two PDE models using the `ClimaCore.jl` backends.

# # Model 1
# Model 1 represents a simplified atmosphere (atm)
# and solves the [heat
# equation](https://en.wikipedia.org/wiki/Heat_equation)
# in a one-column domain:

# ``
# \frac{∂ T}{∂ t} + ∇ ⋅ (-μ ∇T) = 0
# ``

#    with top and bottom boundary conditions set to fixed-temperature (non-zero Dirichlet) and fixed-flux (non-zero Neumann) conditions, respectively:

# ``
#    T_{top}  = 280 K, \,\,\,\,\,\,\, \frac{∂ T_{bottom}}{∂ t} = - ∇ F_{sfc}
# ``

# where
#  - `t` is time
#  - `μ` is the thermal diffusivity
#  - `T` is the temperature
#  - `F_sfc` is the thermal boundary flux (see below for calculation)

# # Model 2
# Model 2 assumes a simplified soil (lnd) domain to be a slab, represented by an ODE for its surface temperature:

# ``
#    \frac{dT_{sfc}}{dt} = - (F_{accumulated} + G ) / h_{lnd}
# ``

#    where

# ``
#    F_{accumulated} = {F_{integrated}} / Δt_{coupler}
# ``

# where
#  - `Δt_coupler` is the duration of the coupling cycle
#  - `T_sfc` is the temperature
#  - `h_lnd` is the slab thickness
#  - `F_integrated` and `F_accumulated` thermal boundary fluxes, respectively (see below for calculation)
#  - `G` represents soil physics, which we assume to be 0 for now

# # Coupling and Flux Calculation

# We use this Model 1 (usually this is done by the model with the shortest timestep)
# to calculate and accumulate the downward surface fluxes, `F_sfc`:

# ``
#    F_{sfc} = - λ (T_{sfc} - T1)
# ``

# ``
#    d(F_{integrated})/dt  = F_{sfc}
# ``

# where
#  - `T1` is the atm temperature near the surface (here assumed equal to the first model level)
#  - `λ` a constant relaxation timescale

# The atmosphere's surface flux, `atm_F_sfc`, is calculated by the atmosphere's timestepping.
# It is then passed to the land model via the coupler. The land model's `lnd_F_sfc` gets updated
# with this value from the atmosphere, but the surface flux is not calculated by the land model itself.

# Note that in a more realistic setup the above equations would be weighted by their domains' densities and thermal heat capacities, so that the thermal flux would have the units of W m$^{-2}$.
# Here we assume these are unity for both domains.

# Sequential coupling has the following steps:
# 1) pre-Model 1: supply Model 1 with `T_sfc` for the `F_sfc` calculation; reset `F_integrated` to zero
# 2) run Model 1: step forward for all Model 1 timesteps within one coupling cycle using `F_sfc` as the bottom boundary condition; accumulate `F_integrated` at each (sub-)step
# 3) post-Model 1: pass `F_integrated` into coupler and convert to `F_accumulated` for the correct units.
# 4) pre-Model 2: supply Model 2 with `F_accumulated`
# 5) run Model 2: step forward for all Model 2 timesteps within one coupling cycle;
# 6) post-Model 2: state variable, `T_sfc` of Model 2 into coupler.
# 7) repeat steps 1-6 for all coupling timesteps.

# # Implementation

# ## [Loading Packages](@id Loading-code-coupler)

# First, we'll load our pre-requisites:
#  - load CliMA packages under development - you may need to add unregistered packages in Pkg, e.g.:
## import Pkg; Pkg.add(url="https://github.com/CliMA/ClimaCore.jl",rev="main")

#  - load external packages:
import SciMLBase

# load CliMA packages:
import ClimaCore as CC
import ClimaTimeSteppers as CTS

# Load utilities for coupling
include("../CoupledSims/coupled_sim.jl")

#src ## Setup Logging Information
const CI = !isnothing(get(ENV, "CI", nothing)) #src

# ## Define Parameters
#  - Global Constants
const FT = Float64;

#  - Experiment-specific Parameters
parameters = (
    ## atmos parameters
    zmin_atm = FT(0.0), # height of atm stack bottom [m]
    zmax_atm = FT(1.0), # height of atm stack top [m]
    n = 15,  # number of elements in atm stack
    μ = FT(0.0001), # diffusion coefficient [m^2 / s]
    T_top = FT(280.0), # fixed temperature at the top of the domain_atm [K]
    T_atm_ini = FT(280.0), # initial condition of at temperature (isothermal) [K]
    ## slab parameters
    h_lnd = FT(0.5), # depth of slab layer [m]
    T_lnd_ini = FT(260.0), # initial condition of at temperature (isothermal) [K]
    ## coupling parameters
    λ = FT(1e-5), # transfer coefficient
)


# ## Define Model Functions

# - Model 1 (atm) Equations
"""
    ∑tendencies_atm!(du, u, (parameters, T_sfc), t)

Heat diffusion equation
    dT/dt =  ∇ μ ∇ T
    where
        T  = 280 K              at z = zmax_atm
        dT/dt = - ∇ F_sfc       at z = zmin_atm

We also use this model to calculate and accumulate the downward surface fluxes, F_sfc:
    F_sfc = - λ * (T_sfc - T1)
    d(F_integrated)/dt  = F_sfc
    where
        F_integrated is reset to 0 at the beginning of each coupling cycle
        T1 = atm temperature near the surface (here assumed equal to the first model level)
"""
function ∑tendencies_atm!(du, u, cache, t)
    parameters, T_sfc = cache
    T_atm = u.T_atm
    F_sfc = calculate_flux(T_sfc[1], parent(T_atm)[1], parameters)

    ## set boundary conditions
    C3 = CC.Geometry.WVector
    # note: F_sfc is converted to a Cartesian vector in direction 3 (vertical)
    bcs_bottom = CC.Operators.SetValue(C3(F_sfc))
    bcs_top = CC.Operators.SetValue(C3(FT(0)))

    ## gradient and divergence operators needed for diffusion in tendency calc.
    ᶠgradᵥ = CC.Operators.GradientC2F()
    ᶜdivᵥ = CC.Operators.DivergenceF2C(bottom = bcs_bottom, top = bcs_top)

    ## perform tendency calculations
    # note: `atm_F_sfc` is a prognostic variable (i.e. is timestepped) in the atmosphere model
    @. du.T_atm = ᶜdivᵥ(parameters.μ * ᶠgradᵥ(T_atm)) # dT/dt
    @. du.atm_F_sfc = -F_sfc[1] # d(F_integrated)/dt
end

# - Model 2 (lnd) Equations
"""
    ∑tendencies_lnd!(dT_sfc, T_sfc, (parameters, F_accumulated), t)

Slab layer equation
    lnd d(T_sfc)/dt = - (F_accumulated + G) / h_lnd
    where
        F_accumulated = F_integrated / Δt_coupler
"""
function ∑tendencies_lnd!(dT_sfc, T_sfc, (parameters, F_accumulated), t)
    G = 0.0 # place holder for soil dynamics
    @. dT_sfc = (-F_accumulated + G) / parameters.h_lnd
end

# - Surface Flux Calculation (coarse bulk formula)
calculate_flux(T_sfc, T1, parameters) = -parameters.λ * (T_sfc - T1);

# - Coupler Communication Functions
# These functions export / import / transform variables
# These functions are now just place holders for coupler transformations (e.g. regridding, masking, etc)
coupler_get_(x) = x;
coupler_put_(x) = x;

# ## Model Initialization
# - initialize atm model domain and grid
domain_atm = CC.Domains.IntervalDomain(
    CC.Geometry.ZPoint{FT}(parameters.zmin_atm),
    CC.Geometry.ZPoint{FT}(parameters.zmax_atm);
    boundary_names = (:bottom, :top),
);
context = CC.ClimaComms.context()
mesh_atm = CC.Meshes.IntervalMesh(domain_atm, nelems = parameters.n); # struct, allocates face boundaries to 5,6: atmos
device = CC.ClimaComms.device(context)
center_space_atm = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh_atm) # collection of the above, discretises space into FD and provides coords

# - initialize prognostic variables, either as ClimaCore's Field objects or as Arrays
T_atm_0 = CC.Fields.ones(FT, center_space_atm) .* parameters.T_atm_ini; # initiates a spatially uniform atm progostic var
T_lnd_0 = [parameters.T_lnd_ini]; # initiates lnd progostic var
ics = (; atm = T_atm_0, lnd = T_lnd_0)

# - specify timestepping information
stepping = (;
    Δt_min = 0.02,
    timerange = (0.0, 6.0),
    Δt_coupler = 0.1,
    odesolver = CTS.ExplicitAlgorithm(CTS.RK4()),
    nsteps_atm = 8, # number of timesteps of atm per coupling cycle
    nsteps_lnd = 1, # number of timesteps of lnd per coupling cycle
);

# ## Define the sequential coupling loop
function coupler_solve!(stepping, ics, parameters)
    ## access start time, end time, and timesteps
    Δt_min = stepping.Δt_min
    Δt_coupler = stepping.Δt_coupler
    t_start = stepping.timerange[1]
    t_end = stepping.timerange[2]

    ## init coupler fields
    coupler_F_sfc = [0.0]
    coupler_T_lnd = copy(ics.lnd)

    ## atmos copies of coupler variables
    atm_T_lnd = copy(coupler_T_lnd)
    atm_F_sfc = copy(coupler_F_sfc)

    ## SETUP ATMOS
    ## put all prognostic variable arrays into a vector and ensure that solve can partition them
    T_atm = ics.atm
    Y_atm = CC.Fields.FieldVector(T_atm = T_atm, atm_F_sfc = atm_F_sfc)

    ode_algo = stepping.odesolver
    ode_function = CTS.ClimaODEFunction(T_exp! = ∑tendencies_atm!)

    # `saveat` determines the frequency with which the solution is saved
    # `adaptive` determines whether the solver should use an adaptive timestep
    prob_atm = SciMLBase.ODEProblem(
        ode_function,
        Y_atm,
        (t_start, t_end),
        (parameters = parameters, T_sfc = atm_T_lnd),
    )
    integ_atm = SciMLBase.init(
        prob_atm,
        ode_algo,
        dt = Δt_min,
        saveat = t_start:(10Δt_min):t_end,
        adaptive = false,
    )

    ## land copies of coupler variables
    # note: `lnd_F_sfc` is a diagnostic variable (i.e. is not timestepped) in the land model
    T_lnd = ics.lnd
    lnd_F_sfc = copy(coupler_F_sfc)

    ## SETUP LAND
    ode_function = CTS.ClimaODEFunction(T_exp! = ∑tendencies_lnd!)

    prob_lnd =
        SciMLBase.ODEProblem(ode_function, T_lnd, (t_start, t_end), (parameters, lnd_F_sfc))
    integ_lnd = SciMLBase.init(
        prob_lnd,
        ode_algo,
        dt = Δt_min,
        saveat = t_start:(10Δt_min):t_end,
        adaptive = false,
    )


    ## coupler stepping
    for t in ((t_start + Δt_coupler):Δt_coupler:t_end)

        ## STEP ATMOS
        ## pre_atmos
        integ_atm.p[2] .= coupler_get_(coupler_T_lnd) # integ_atm.p is the parameter vector of an ODEProblem from SciMLBase
        integ_atm.u.atm_F_sfc .= 0.0 # surface flux to be accumulated

        ## run atmos
        ## NOTE: use (t - integ_atm.t) here instead of Δt_coupler to avoid accumulating roundoff error in our timestepping.
        SciMLBase.step!(integ_atm, t - integ_atm.t, true)

        ## post_atmos
        ## retrieve the updated surface flux from the atmosphere model
        coupler_F_sfc .= coupler_put_(integ_atm.u.atm_F_sfc) / Δt_coupler

        ## STEP LAND
        ## pre_land
        ## provide the updated surface flux to the land model
        lnd_F_sfc .= coupler_get_(coupler_F_sfc)

        ## run land
        SciMLBase.step!(integ_lnd, t - integ_lnd.t, true)

        ## post land
        coupler_T_lnd .= coupler_put_(integ_lnd.u) # update T_sfc
    end

    return integ_atm, integ_lnd
end;

# ## Run the Coupler Model Simulation
integ_atm, integ_lnd = coupler_solve!(stepping, ics, parameters);
sol_atm, sol_lnd = integ_atm.sol, integ_lnd.sol;

# ## Postprocessing and Visualization

# Each integrator output (`sol_atm`, `sol_lnd`), contains the SciMLBase variable `.u` (the name is hard coded).
# If multiple prognostic variables were solved for, `u` will contain all of those.
# `parent()` accesses the `Field` values.
# So, for example, the structure of `u` from the atmosphere model is:
# `parent(sol_atm.u[<time-index>].T_atm)[<z-index>,<variable-index>]`
# and
# `parent(sol_atm.u[<time-index>].atm_F_sfc)[<z-index>,<variable-index>]`

ENV["GKSwstype"] = "nul"
import Makie
import CairoMakie

ARTIFACTS_DIR = joinpath("experiments/ClimaCore/output/heat-diffusion/artifacts")
mkpath(ARTIFACTS_DIR)

# - Vertical profile at start and end
t0_ = parent(sol_atm.u[1].T_atm)[:, 1];
tend_ = parent(sol_atm.u[end].T_atm)[:, 1];
z_centers = parent(CC.Fields.coordinate_field(center_space_atm))[:, 1];
profile_fig = Makie.Figure();
profile_ax = Makie.Axis(
    profile_fig[1, 1],
    xlabel = "T (K)",
    ylabel = "z (m)",
    title = "Atmos profile at start & end of Simulation",
)
Makie.lines!(profile_ax, t0_, z_centers, label = "t=0"; linewidth = 2)
Makie.lines!(profile_ax, tend_, z_centers, label = "t=end"; linewidth = 2)
Makie.axislegend(profile_ax, position = :lt)
Makie.save(joinpath(ARTIFACTS_DIR, "atmos_profile.png"), profile_fig)

# - Conservation: absolute "energy" of both models with time
# convert to the same units (analogous to energy conservation, assuming that is both domains density=1 and thermal capacity=1)
lnd_sfc_u_t = [u[1] for u in sol_lnd.u] .* parameters.h_lnd;
atm_sum_u_t =
    [sum(parent(u.T_atm)[:]) for u in sol_atm.u] .* (parameters.zmax_atm - parameters.zmin_atm) ./
    parameters.n;
v1 = lnd_sfc_u_t .- lnd_sfc_u_t[1];
v2 = atm_sum_u_t .- atm_sum_u_t[1];
component_energy_fig = Makie.Figure();
component_energy_ax = Makie.Axis(
    component_energy_fig[1, 1],
    xlabel = "time (s)",
    ylabel = "energy flux (J / m2)",
    title = "Component Model Energy during Simulation",
)
Makie.lines!(component_energy_ax, sol_lnd.t, v1, label = "lnd"; linewidth = 2)
Makie.lines!(component_energy_ax, sol_lnd.t, v2, label = "atm"; linewidth = 2)
Makie.lines!(component_energy_ax, sol_lnd.t, v1 .+ v2, label = "tot"; linewidth = 2)
Makie.axislegend(component_energy_ax, position = :lt)
Makie.save(joinpath(ARTIFACTS_DIR, "component_energy.png"), component_energy_fig)

# - Conservation: relative error with time
mean(arr) = sum(arr) / length(arr)
total = atm_sum_u_t + lnd_sfc_u_t;
rel_error = abs.(total .- total[1]) / mean(total);
energy_conservation_fig = Makie.Figure();
energy_conservation_ax = Makie.Axis(
    energy_conservation_fig[1, 1],
    xlabel = "simulation time (s)",
    ylabel = "relative error",
    title = "Total Energy Conservation",
)
Makie.lines!(energy_conservation_ax, sol_lnd.t, rel_error, linewidth = 2)
Makie.save(joinpath(ARTIFACTS_DIR, "energy_conservation.png"), energy_conservation_fig)

#src # - Animation
#src anim = Plots.@animate for u in sol_atm.u
#src     Plots.plot(u.T_atm, xlim=(220,280))
#src end
#src Plots.mp4(anim, joinpath(path, "heat.mp4"), fps = 10)
#src
#src function linkfig(figpath, alt = "")
#src     # buildkite-agent upload figpath
#src     # link figure in logs if we are running on CI
#src     if get(ENV, "BUILDKITE", "") == "true"
#src         artifact_url = "artifact://$figpath"
#src         print("\033]1338;url='$(artifact_url)';alt='$(alt)'\a\n")
#src     end
#src end
#src
#src linkfig(ARTIFACTS_DIR * "/heat_end.png", "Heat End Simulation")
