# Atmos Single Column & Slab Land


            |============================== z = 1e4 m  
        A   |
        T   |              ...
        M   |
        O   |------------------------------ z = 1e3 m    
        S   |
            |
            |============================== z = 0 m    


       LAND |============================== z = 0 m  


## **Atmos Component**

### 1. Formulation
We're using the ClimateMachine.jl's Atmos SingleStack configuration. In this simplified case, we are solving the 5 diffusion equations: 1 for density, 3 for the velocity components, and 1 for a tracer. 

$$
\frac{\partial ρ}{\partial t} = - ∇ \cdot (\vec{u} ρ + μ ∇ ρ) \qquad 
\frac{\partial \vec{u}}{\partial t} = - ∇ \cdot (\vec{u} \vec{u} + ν ∇ \vec{u}) \qquad 
\frac{\partial θ}{\partial t} = - ∇ \cdot (\vec{u} \theta + κ ∇ θ) \qquad 
$$
with μ, ν and κ representing the respective viscosity and diffusivity coefficients. The boundary conditions are:
- Upper (external) boundary:
    - ρ: `Impenetrable()` 
    - u: `FreeSlip()`
    - θ: `Insulating()`
- Lower (coupled) boundary:
    - ρ: `Impenetrable()`  
    - u: `FreeSlip()`
    - θ: `CoupledPrimary()`: so that at the boundary $\frac{\partial θ}{\partial t} = ... - ∇ \cdot F_{tot}$ where $F_{tot}$ is the surface flux calculated by `calculate_land_sfc_fluxes(_...)` using Atmos and Land properties (see below) 


    NB: See the [boundary_conditions_doc.md](../../docs/src/boundary_conditions_doc.md) for more details on BCs and their implementation.

### 2. Implementation notes
- We followed the Bickley jet the ClimateMachine.jl's RQA interface, except the boundary condition for ρθ. In this case, we impose a non-zero normal diffusive flux, which is added to the total normal flux, `fluxᵀn` (also referred to as `local_flux` in the code):

        function numerical_boundary_flux_second_order!(
            numerical_flux::Union{PenaltyNumFluxDiffusive},
            bctype::CoupledPrimaryBoundary,
            balance_law::ModelSetup,
            fluxᵀn::Vars, _...)
            F_tot = calculate_land_sfc_fluxes(balance_law, state_prognostic⁻, state_auxiliary⁻, t)  # W/m^2
            fluxᵀn.ρθ += - F_tot  / balance_law.parameters.cp_d

        end

## **Slab Land Component**

For prototyping the coupler, we parameterize the land as a 0-dimensional slab covering the whole globe, analogous to the slab ocean with no water. 

### 1. Formulation
The evolution of the prognostic variable, surface soil temperature ($T_{sfc}$), can be written in the units of heat flux (W m$^{-2}$):

$$
\rho_s h_s c_s \frac{\partial T_{sfc}}{\partial t} = F_{tot} + G = R_{SW} - R_{LW} - SH - LH + G
$$

  
where

- Sensible heat flux at surface
    $$SH = cp_d g_a (T_{sfc} - T_a) = cp_d g_a (T_{sfc} - T_a)$$
- Latent heat flux at surface
    $$LH = \lambda g_w (q_{sat}(T_{sfc},p_{sfc}) - q_a)$$
- Downward-pointing short-wave radiative flux
    $$R_{SW} = F_{sfc}^d - F_{sfc}^u = (1-\alpha)\tau F_{sol} * (1 + \sin(2π τ_d^{-1} t))$$
- Upward-pointing net long-wave flux
    $$R_{LW} = F_{sfc}^u - F_{sfc}^d = \epsilon \sigma T_{sfc}^4 - \epsilon F_{a}$$
- Soil storage
    $$G = \kappa_s (T_{sfc}-T_h)/h_s$$

For description and values of the parameters above, please see [parameters_initialconditions.jl](parameters_initialconditions.jl). 

### 2. Implementation
- For the slab land, we use explicit time stepping:
$$T_{sfc}^{n+1} = T_{sfc}^n +\frac{(F_{tot}^{n} + G^{n})}{\rho_s c_s h_s}\Delta t$$
- This is implemented as:

        @inline function source!(
                model::SlabLandModelSetup, source::Vars, state::Vars, gradflux::Vars, aux::Vars, _...)
            p = model.parameters
            G    = p.κ_s * (state.T_sfc - p.T_h) / p.h_s # simple soil physics
            source.T_sfc = - (aux.F_ρθ_prescribed + G) / (p.ρ_s * p.c_s * p.h_s)
        end

## **Coupler**

### 1. Formulation
- we're using our sequential coupler, which is described [here](https://clima.github.io/CouplerMachine/dev/timestepping/)
### 2. Surface Flux Handling
1. Flux calculation
    - function `calculate_land_sfc_fluxes(_...)` is defined in `coupler/surface_fluxes.jl`, and is used for:
        - the Atmos boundary condition calculation (see above)
        - for surface flux accumulation (next point). 
    - it is called by Atmos, which is assumed to be the model with the shortest timestep. This means that the fluxes are only calculated once and are consistent with Atmos's time stepping.

2. Flux accumulation
    - In order to integrate fluxes in time consistently with the Atmos time stepper, we define a new prognostic variable `F_ρθ_accum` whose source function accumulates the fluxes by integrating them in time, so that $F\_ρθ\_accum = \int F_{tot} dt$.

            @inline function source!(atmos_model::ModelSetup, _...)
                source.F_ρθ_accum = calculate_land_sfc_fluxes(atmos_model, _...) 
            end

3. Flux / state storage and communication
    - the coupler stores fields:
        - `BoundaryEnergyFlux` which collects `F_ρθ_accum / Δt` calculated and accumulated in Atmos during the Atmos part of the coupling cycle (with $Δt$ denoting the period of the coupling cycle), and is called by land for use in its $T_{sfc}$ equation
        - `LandSurfaceTemperature` which collects the state of the slab Land model, namely $T_{sfc}$, and is called by atmos when calculating surface fluxes


## **Tests**
- Dry heat diffusion
    - Simplifications
        - $g_a$ is constant
        - (μ, ν, $\kappa$) = (0, 0, 1e-5)
        - $\vec{u} = (0,0,0)$
        - $G = 0$ and $F_a = 0$
    - [results](https://docs.google.com/document/d/1JKK8wFKPq3Jo3D4flXZiY3WAUPqYE19pXRvwJgpDwjw/edit)
- Dry advection-diffusion
    - Simplifications
        - $g_a$ is constant
        - (μ, ν, $\kappa$) = (0, 0, 1e-5)
        - $\vec{u} = (1,0,0)$
        - $G \neq 0$ and $F_a = 0$
    - [results](https://docs.google.com/document/d/1JKK8wFKPq3Jo3D4flXZiY3WAUPqYE19pXRvwJgpDwjw/edit)

## **Pipeline**

For this simple implementation, we assume:
-  a dry setup 
- atmosphere is only governed by the advection-diffusion
- aerodynamic conductance (i.e. $g_a = |u| c_D$, with $c_D$ representing the drag coefficient) is constant. 

Using explicit time stepping the formulation is implemented as follows:

1) Both models' states are initiallised according to the user-defined ICs 
2) Land model performs its first part of the coupling cycle, with zero surface fluxes and saves its state into the coupler field `LandSurfaceTemperature`
3) Coupler transforms its field `LandSurfaceTemperature` into Atmos field `auxiliary.T_sfc`
4) Atmos performs all its timesteps within the coupling cycle, with the boundary conditions at the coupled boundary contributing to the total normal flux with the surface flux coming from the land (`fluxᵀn.ρθ += - F_tot / cp_d`, so that $\partial_t \rho \theta = ... + \nabla_z \sdot F_{tot}/cp_d$)
5) The same flux `F_tot` is calculated in the `source!` function and saved as `state.F_ρθ_accum`, which ensures that the flux will be integrated and accumulated in time, consistent with the Atmos time stepping. 
6) Coupler converts Atmos field `state.F_ρθ_accum` to the coupler field `BoundaryEnergyFlux` 
7) Coupler sends its field `BoundaryEnergyFlux` and transforms it into Land field `auxiliary.F_ρθ_prescribed`
8) Land performs its own initialization and physics (`G`), which is then added to the  `auxiliary.F_ρθ_prescribed` (corresponding to $F_{tot}$ above) in the `source!` function to solve:
$$T_s^{n+1} = T_s^n +\frac{F_{tot}^{n...}}{\rho_s c_s h_s}\Delta t$$
The $^{n}$ fields correspond to the previous Land timestep.
9) Coupler converts Land field `state.T_sfc` to the coupler field `LandSurfaceTemperature` 
10) Atmos `state.F_ρθ_accum` is reset to 0 and the cycle repeats.

## References:
- [Bonan 2019 book](https://www.cambridge.org/us/academic/subjects/earth-and-environmental-science/climatology-and-climate-change/climate-change-and-terrestrial-ecosystem-modeling?format=HB&isbn=9781107043787)
- https://www.cesm.ucar.edu/models/atm-cam/docs/description/node29.html
- http://www.met.reading.ac.uk/~swrhgnrj/teaching/MT23E/mt23e_notes.pdf