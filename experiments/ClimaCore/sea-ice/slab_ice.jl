# # Heat Equation + Slab Tutorial

#  - load external packages:
import LinearAlgebra, UnPack
import ClimaCore: Fields, Domains, Topologies, Meshes, DataLayouts, Operators, Geometry, Spaces

push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", ".."))
using ClimaCoupler

using Base: show_supertypes
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33

using Logging: global_logger
using TerminalLoggers: TerminalLogger

using RecursiveArrayTools

using OrdinaryDiffEq

using Statistics

#src ## Setup Logging Information
global_logger(TerminalLogger())
const CI = !isnothing(get(ENV, "CI", nothing))

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
    ## slab ice parameters
    ρc_ml = FT(4e6), # density times heat transfer coefficient for mixed layer [J / m2 / K ]
    F0_base = FT(120), # ice base transfer coefficient [W / m2 / K]
    T_base = FT(273.16), # ice base temperature [K] 
    L_ice = FT(3e8), # latent heat coefficient for ice [J / m3]
    h_ml = FT(1), # mixed layer depth [m] 
    T_freeze = FT(273.16), # temperature at freezing point [K] 
    k_ice = FT(2), # thermal conductivity of ice [W / m / K]
    ## coupling parameters 
    λ = FT(0.001), # coupling coefficient (J / K / s / m2)
)

# ## Define Model Functions

# - Model 1 (atm) Equations
"""
    ∑tendencies_atm!(du, u, (parameters, T_sfc), t)

Heat diffusion equation
    dT/dt =  ∇ μ ∇ T
    where
        (- μ ∇ T) = 0           at z = zmax_atm
        (- μ ∇ T) = F_sfc       at z = zmin_atm

We also use this model to calculate and accumulate the downward surface fluxes, F_sfc:
    F_sfc = - λ * (T_sfc - T1) 
    d(F_integrated)/dt  = F_sfc
    where
        F_integrated is reset to 0 at the beginning of each coupling cycle
        T1 = atm temperature near the surface (here assumed equal to the first model level)
"""
function ∑tendencies_atm!(du, u, p, t)

    FT = eltype(u)
    T = u.T # u.x = vector of prognostic variables from DifferentialEquations
    F_sfc = calculate_flux(p.Ya.T_sfc[1], parent(T)[1], p.p)

    ## set BCs

    bcs_bottom = Operators.SetValue(Geometry.WVector(F_sfc ./ p.p.λ))
    bcs_top = Operators.SetValue(Geometry.WVector(FT(0)))

    gradc2f = Operators.GradientC2F()
    gradf2c = Operators.DivergenceF2C(bottom = bcs_bottom, top = bcs_top) # # Neumann BC (face-to-center)
    ## tendency calculations
    @. du.T = gradf2c(p.p.μ * gradc2f(T)) # dT/dt
    du.F .= -F_sfc[1] # d(F_integrated)/dt

end

get_∂F_atm∂T_sfc(p) = -p.λ # TODO: in GCM this will be: 4 σ T_sfc^3 + ||uh|| ρ cp (C_H + L(q_sat1 - q_sat) * dT^(-1))

# - Model 2 (ice) Equations
"""
    solve_ice!(dT_sfc, T_sfc, (parameters, F_accumulated), t)

slab RHS with an implicit solve ice and explicit (forward Euler) solve for ocean

"""
function solve_ice!(integ, Δt)

    Y = integ.u
    Ya = integ.p.Ya
    p = integ.p.p

    FT = eltype(Y)
    ocean_qflux = FT(0)

    # prognostic
    T_ml = Y.T_ml
    T_sfc = Y.T_sfc
    h_ice = Y.h_ice

    # auxiliary
    F_atm = Ya.F_atm

    ∂F_atm∂T_sfc = get_∂F_atm∂T_sfc(p) # this will be passed from atmos/SF.jl

    ΔT_ml = similar(F_atm)
    Δh_ice = similar(F_atm)
    F_ice = similar(F_atm)
    ΔT_sfc = similar(F_atm)

    # ice thickness and mixed layer temperature changes due to atmosphereic and ocean fluxes 
    if h_ice[1] > 0 # ice-covered 
        @. ΔT_ml = -(p.F0_base * (T_ml - p.T_base) + ocean_qflux) * Δt / (p.h_ml * p.ρc_ml)
        @. Δh_ice = (F_atm - p.F0_base * (T_ml - p.T_base)) * Δt / p.L_ice
    else # ice-free
        @. ΔT_ml = -(F_atm + ocean_qflux) * Δt / (p.h_ml * p.ρc_ml)
        @. Δh_ice = 0
    end

    # adjust if transition to ice-covered 
    if (T_ml[1] + ΔT_ml[1] < p.T_freeze)
        @. Δh_ice = Δh_ice - (T_ml + ΔT_ml - p.T_freeze) * (p.h_ml * p.ρc_ml) / p.L_ice
        @. ΔT_ml = p.T_freeze - T_ml
    end

    # adjust if transition to ice-free
    if ((h_ice[1] > 0) & (h_ice[1] + Δh_ice[1] <= 0))
        @. ΔT_ml = ΔT_ml - (h_ice + Δh_ice) * p.L_ice / (p.h_ml * p.ρc_ml)
        @. Δh_ice = -h_ice
    end

    # solve for T_sfc
    if (h_ice[1] + Δh_ice[1] > 0) #  surface is ice-covered
        # if ice covered, solve implicity (for now one Newton iteration: ΔT_s = - F(T_s) / dF(T_s)/dT_s )
        @. F_ice = p.k_ice / (h_ice + Δh_ice) * (p.T_base - T_sfc)
        @. ΔT_sfc = (-F_atm + F_ice) / (p.k_ice / (h_ice + Δh_ice) + ∂F_atm∂T_sfc)
        if (T_sfc[1] + ΔT_sfc[1] > p.T_freeze)
            @. ΔT_sfc = p.T_freeze - T_sfc
        end
        # surface is ice-covered, so update T_sfc as ice surface temperature
        @. Y.T_sfc += ΔT_sfc
    else # ice-free, so update T_sfc as mixed layer temperature
        @. Y.T_sfc = T_ml + ΔT_ml
    end

    # update state
    @. Y.T_ml += ΔT_ml
    @. Y.h_ice += Δh_ice

    @. Ya.ice_mask = h_ice[1] > 0 ? h_ice[1] * FT(1) : h_ice[1] * FT(0)

    return
end

function ∑tendencies_ice_stub(du, u, p, t)
    dY = du
    FT = eltype(dY)

    solve_ice!((; u = u, p = p), p.Δt) # timestepping outside of DeffEq (but DeffEq still used here for saving vars in `integ.sol`)

    @. dY.T_ml = FT(0)
    @. dY.h_ice = FT(0)
    @. dY.T_sfc = FT(0)

end

# - Surface Flux Calculation (coarse bulk formula)
calculate_flux(T_sfc, T1, parameters) = -parameters.λ * (T_sfc - T1);

# ## Model Initialization
# - initialize atm model domain and grid
domain_atm = Domains.IntervalDomain(
    Geometry.ZPoint{FT}(parameters.zmin_atm),
    Geometry.ZPoint{FT}(parameters.zmax_atm);
    boundary_tags = (:bottom, :top),
);
mesh_atm = Meshes.IntervalMesh(domain_atm, nelems = parameters.n); # struct, allocates face boundaries to 5,6: atmos
center_space_atm = Spaces.CenterFiniteDifferenceSpace(mesh_atm); # collection of the above, discretises space into FD and provides coords

# - initialize prognostic variables, either as ClimaCore's Field objects or as Arrays
T_atm_0 = Fields.ones(FT, center_space_atm) .* FT(265);
Y_atm = Fields.FieldVector(T = T_atm_0, F = [FT(0)])
Y_ice = Fields.FieldVector(T_sfc = [parameters.T_freeze], h_ice = [FT(0)], T_ml = [parameters.T_freeze])

# - initialize auxiliary variables
ics_aux = (; atm = (; T_sfc = copy(Y_ice.T_sfc)), ice = (; ice_mask = [FT(0)], F_atm = [FT(0)]))

# - specify timestepping information
stepping = (;
    Δt_min = 0.01,
    timerange = (0.0, 10.0),
    Δt_coupler = 0.01,
    odesolver = Euler(), #SSPRK33(),
    nsteps_atm = 1, # number of timesteps of atm per coupling cycle
    nsteps_ice = 1, # number of timesteps of ice per coupling cycle
);

# ## Define the sequential coupling loop
function coupler_solve!(stepping, ics_aux, parameters)
    t = 0.0
    Δt_min = stepping.Δt_min
    Δt_coupler = stepping.Δt_coupler
    t_start = stepping.timerange[1]
    t_end = stepping.timerange[2]

    ## SETUP ATMOS
    ## put all prognostic variable arrays into a vector and ensure that solve can partition them
    prob_atm = ODEProblem(∑tendencies_atm!, Y_atm, (t_start, t_end), (; p = parameters, Ya = ics_aux.atm))
    integ_atm = init(prob_atm, stepping.odesolver, dt = Δt_min, saveat = 1 * Δt_min)

    ## SETUP ICE
    prob_ice =
        ODEProblem(∑tendencies_ice_stub, Y_ice, (t_start, t_end), (; p = parameters, Ya = ics_aux.ice, Δt = Δt_coupler))
    integ_ice = init(prob_ice, Euler(), dt = Δt_coupler, saveat = 1 * Δt_coupler)

    ## SETUP COUPLER
    coupler = CouplerState()
    # pointers to the coupled fields are passed to the coupler:
    coupler_add_field!(coupler, :T_sfc_ice, integ_ice.u.T_sfc)
    coupler_add_field!(coupler, :F_atm, integ_atm.u.F)

    ## coupler stepping
    for t in (t_start:Δt_coupler:t_end)

        ## STEP ATMOS
        ## atmos pull!
        integ_atm.u.F .= [0.0] # reset surface flux to be accumulated; need to do this AFTER ice pulls from coupler
        integ_atm.p.Ya.T_sfc .= coupler_get(coupler, :T_sfc_ice) # integ_atm.p is the parameter vector of an ODEProblem from DifferentialEquations

        ## run atmos
        ## NOTE: use (t - integ_atm.t) here instead of Δt_coupler to avoid accumulating roundoff error in our timestepping.
        OrdinaryDiffEq.step!(integ_atm, Δt_coupler, true)

        # ## atmos push! (not strictly necessary since this is a pointer)
        coupler_put!(coupler, :F_atm, integ_atm.u.F)

        ## STEP ICE
        ## ice pull!
        integ_ice.p.Ya.F_atm .= coupler_get(coupler, :F_atm) / Δt_coupler
        Δt_coupler_ice = Δt_coupler

        ## run ice
        OrdinaryDiffEq.step!(integ_ice, Δt_coupler_ice, true)

        # ## ice push!
        coupler_put!(coupler, :T_sfc_ice, integ_ice.u.T_sfc)

    end

    return integ_atm, integ_ice
end;

# ## Run the Coupler Model Simulation
integ_atm, integ_ice = coupler_solve!(stepping, ics_aux, parameters);
sol_atm, sol_ice = integ_atm.sol, integ_ice.sol;

# ## Postprocessing and Visualization

# Each integrator output (`sol_atm`, `sol_ice`), contains the DifferentialEquations variable `.u` (the name is hard coded).

ENV["GKSwstype"] = "nul"
import Plots
Plots.GRBackend()

show_plots = isdefined(Main, :SHOWPLOTS) ? SHOWPLOTS : true

path = joinpath(dirname(@__FILE__), "images/")
mkpath(path);

# - Vertical profile at start and end
t0_ = parent(sol_atm.u[1].T)[:, 1];
tend_ = parent(sol_atm.u[end].T)[:, 1];
z_centers = parent(Fields.coordinate_field(center_space_atm))[:, 1];
show_plots ?
Plots.png(
    Plots.plot(
        [t0_ tend_],
        z_centers,
        title = "model 1: atm",
        labels = ["t=0" "t=end"],
        xlabel = "T (K)",
        ylabel = "z (m)",
    ),
    joinpath(path, "tc1_f1.png"),
) : nothing
# ![](images/tc1_f1.png)

# - Conservation: absolute "energy" of both models with time
# convert to the same units (analogous to energy conservation, assuming that is both domains density=1 and thermal capacity=1)
ice_sfc_h_t = [sum(parent(u.h_ice)[:]) for u in sol_ice.u] .* parameters.h_ml * parameters.L_ice;
ice_sfc_u_t = [sum(parent(u.T_sfc)[:]) for u in sol_ice.u] .* parameters.h_ml * parameters.ρc_ml .- ice_sfc_h_t;
atm_sum_u_t =
    [sum(parent(u.T)[:]) for u in sol_atm.u] .* (parameters.zmax_atm .- parameters.zmin_atm) ./ parameters.n .* parameters.λ;
v1 = ice_sfc_u_t .- ice_sfc_u_t[1];
v2 = atm_sum_u_t .- atm_sum_u_t[1];
show_plots ?
Plots.png(
    Plots.plot(
        sol_ice.t,
        [v1 v2 v1 + v2],
        labels = ["ice" "atm" "tot"],
        xlabel = "time (s)",
        ylabel = "pseudo-energy (J / m2)",
    ),
    joinpath(path, "tc1_f2.png"),
) : nothing
# ![](images/tc1_f2.png)

# - Conservation: relative error with time
total = atm_sum_u_t + ice_sfc_u_t;
rel_error = (total .- total[1]) / mean(total);
show_plots ?
Plots.png(
    Plots.plot(sol_ice.t, rel_error, labels = ["tot"], xlabel = "time (s)", ylabel = "relative error"),
    joinpath(path, "tc1_f3.png"),
) : nothing
# ![](images/tc1_f3.png)

#src # - Animation
#src anim = Plots.@animate for u in sol_atm.u
#src     Plots.plot(u.x[1], xlim=(220,280))
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
#src dirname = "heat"
#src linkfig("output/$(dirname)/heat_end.png", "Heat End Simulation")
