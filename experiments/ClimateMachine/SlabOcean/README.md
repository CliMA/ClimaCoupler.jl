# Moist Atmos GCM & Slab Ocean


    Two spherical shells:

                |==z = 3e4 m =====
            A   |                  \\ 
            T   |      ...           
            M   |
            O   |--z = 1e3 m ---  
            S   |                \
                |
                |==z = 0 m ==  
                              \\

           LAND |==z = 0 m ==
                              \\

## **Atmos Component**

We use the ClimateMachine.jl's Atmos GCM moist Held-Suarez configuration, with the RQA interface (and the ESDGModel backend), loosely following [Thatcher & Jablonowski (2016)](https://gmd.copernicus.org/articles/9/1263/2016/gmd-9-1263-2016.pdf). We solve six prognostic equations for $ρ$, $ρ\vec{u}$ (3 components), $ρe$ and $ρq$. We solve the linear part (with respect to ρ and ρu) of this set implicitly in order to alleviate problems arising from acoustic and gravity waves, and the rest is solved explicitly using Runge Kutta timestepping. With this ImEx timestepping, the two equation sets we solve are:

### **1. Full model**
- Equations:

    $$
    \frac{\partial ρ}{\partial t} + ∇ \cdot (\vec{u} ρ ) = - \max(0, ρq - ρ qᵥ) / τ \qquad 
    $$
    $$
    \frac{\partial \vec{ρu}}{\partial t} + ∇ \cdot (\vec{ρu} \otimes \vec{u} + p I ) = - 2 \Omega  \times ρ\vec{v} - ρ \nabla Φ - k_v ρ\vec{v} \qquad  
    $$
    $$
    \frac{\partial ρe}{\partial t} + ∇ \cdot (\vec{u} (ρe + p)) = - (c_{vl}(T - T_0) + Φ) \max(0, ρq - ρ qᵥ) / τ - k_T ρ c_{vd} (T - T_{equil}) \qquad 
    $$
    $$
    \frac{\partial ρq}{\partial t} + ∇ \cdot (\vec{u} ρq ) = - \max(0, ρq - ρ qᵥ)  / τ \qquad 
    $$
    where $I$ is the rank-3 identity matrix, $\vec{v}$ is the horizontal velocity vector, and $T_0 = 273.16 K$. For a more detailed description and values of the parameters above, please see [parameters_initialconditions.jl](parameters_initialconditions.jl). 

- Boundary conditions
    - Upper (external) boundary:
        - ρ: `Impenetrable()` 
        - ρu: `FreeSlip()`
        - ρe, ρq: `Insulating()`
    - Lower (coupled) boundary:
        - ρ: `Impenetrable()`  
        - ρu: `FreeSlip()`
        - ρe, ρq: `CoupledPrimary()`: $\frac{\partial θ}{\partial t} = ... - ∇ \cdot F_{tot}$ where $F_{tot}$ is the surface flux calculated by `calculate_land_sfc_fluxes(_...)` using Atmos and Ocean properties (see below) 

                        function numerical_boundary_flux_second_order!(
                            numerical_flux::Union{PenaltyNumFluxDiffusive},
                            bctype::CoupledPrimaryBoundary,
                            balance_law::ModelSetup,
                            fluxᵀn::Vars, _...)
                            F_tot = calculate_land_sfc_fluxes(balance_law, state_prognostic⁻, state_auxiliary⁻, t)  # W/m^2
                            fluxᵀn.ρθ = - F_tot  / balance_law.parameters.cp_d

                        end

        NB: See the [boundary_conditions_doc.md](../../docs/src/boundary_conditions_doc.md) for more details on BCs and their implementation.

### **2. Linear model**
- Equations

    $$
    \frac{\partial ρ}{\partial t} + ∇ \cdot (\vec{u} ρ) = 0 \qquad 
    \frac{\partial \vec{ρu}}{\partial t} + ∇ \cdot p_L  = - ρ \nabla Φ \qquad 
    \frac{\partial ρe}{\partial t} + ∇ \cdot \frac{(ρeᵣ + pᵣ)  ρu }{ρᵣ} = 0  
    $$

    where the linear pressure is calculated as

    $$
    p_L = (c\_p_d / c\_v_d - 1)  (ρe - ρ  Φ - ρe_{latent} + ρ c\_v_d T_0)
    $$
    $$
    ρe_{latent} = (ρq_{tot} - ρq_{liq})  e\_int_{v0} - ρq_{ice}  (e\_int_{v0} + e\_int_{i0})
    $$

    NB: For a more detailed description and values of the parameters above, please see [parameters_initialconditions.jl](parameters_initialconditions.jl). 
    

-  Boundary conditions
    - ρ: `Impenetrable()`; ρu: `FreeSlip()`; ρe, ρq: `Insulating()` for both boundaries


## **Slab Ocean Component**

For prototyping of the coupler, we parameterize the land as a 0-dimensional slab covering the whole globe. 

### 1. Formulation


$$
\rho_o h_o c_o \frac{\partial T_{sfc}}{\partial t} = F_{tot} + Q = R_{SW} - R_{LW} - SH - LH + Q
$$

  
where

- Sensible heat flux at surface
    $$SH = ρ C_e |u| c_{vd} (T_{sfc} - T_a)$$
- Latent heat flux at surface
    $$LH = ρ C_l |u| L_{Hv0} (q_{sat}(T_{sfc},p_{sfc}) - q_a)$$
- Downward-pointing short-wave radiative flux
    $$R_{SW} = F_{sfc}^d - F_{sfc}^u = (1-\alpha)\tau F_{sol} (1 + \sin(2π τ_d^{-1} t))$$
- Upward-pointing net long-wave flux
    $$R_{LW} = F_{sfc}^u - F_{sfc}^d = \epsilon \sigma \theta_{sfc}^4 - \epsilon F_{a}$$
- Ocean mixing processes are parameterized using a time invariant Q-flux (formulation follows the FMS idealized model)
    $$Q = Q_0\frac{(1-2 \phi^2/L_w^2)}{\cos\phi}exp^{- (\phi^2/L_w^2)}$$

For description and values of the parameters above, please see [parameters_initialconditions.jl](parameters_initialconditions.jl). 

### 2. Implementation
- for the slab ocean, we use explicit time stepping:
$$T_{sfc}^{n+1} = T_{sfc}^n +\frac{(F_{tot}^{n} + Q^{n})}{\rho_o c_o h_o}\Delta t$$
- this is implemented as

        @inline function source!(
                model::SlabOceanModelSetup, source::Vars, state::Vars, gradflux::Vars, aux::Vars, _...)
            p = model.parameters
            Q    =  idealized_qflux(p)
            source.T_sfc = (aux.F_ρθ_prescribed + Q) / (p.ρ_s * p.c_s * p.h_s)
        end

## **Coupler**

### 1. Formulation
- we're using our sequential coupler, which is described [here](https://clima.github.io/ClimaCoupler.jl/dev/timestepping/)
### 2. Surface Flux Handling
1. Flux calculation
    - function `calculate_ocean_sfc_fluxes(_...)` is defined in `coupler/surface_fluxes.jl`, and is used for:
        - the Atmos boundary condition calculation (see above)
        - for surface flux accumulation (next point). 
    - it is called by Atmos, which is assumed to be the model with the shortest timestep. This means that the fluxes are only calculated once and are consistent with Atmos's time stepping.

2. Flux accumulation
    - In order to integrate fluxes in time consistently with the Atmos time stepper, we define a new prognostic variable `F_ρθ_accum` whose source function accumulates the fluxes by integrating them in time, so that $F\_ρθ\_accum = \int F_{tot} dt$.

            @inline function source!(atmos_model::ModelSetup, _...)
                source.F_ρθ_accum = calculate_ocean_sfc_fluxes(atmos_model, _...) 
            end

3. Flux / state storage and communication
    - the coupler stores fields:
        - `BoundaryEnergyFlux` which collects `F_ρθ_accum / Δt` calculated and accumulated in Atmos during the Atmos part of the coupling cycle (with $Δt$ denoting the period of the coupling cycle), and is called by land for use in its $T_{sfc}$ equation
        - `SeaSurfaceTemperature` which collects the state of the slab Land model, namely $T_{sfc}$, and is called by atmos when calculating surface fluxes

## **Tests**
- conservation checks: dry SlabOcean w/o Held-Suarez forcing
    - [results (June 23 section)](https://docs.google.com/document/d/1JKK8wFKPq3Jo3D4flXZiY3WAUPqYE19pXRvwJgpDwjw/edit)
- physical checks: moist SlabOcean with Held-Suarez forcing
    - [results (July 9 section)](https://docs.google.com/document/d/1JKK8wFKPq3Jo3D4flXZiY3WAUPqYE19pXRvwJgpDwjw/edit)

## **Other Notes**
- diffusion can be added once available in the ESDGModel kernels
- currently, there is a bug in LMARS numerical fluxes, leading to conservation breaking at the boundary; as an alternative, Roe fluxes can be used (but the moist version is yet to be implemented, and in the dry case Roe fluxes have been found more unstable due to their less extreme smoothing at boundaries)  

## **References**
- [CESM slab ocean model](https://www.cesm.ucar.edu/models/atm-cam/docs/description/node29.html)
- [GFDL mixed-layer model](https://www.gfdl.noaa.gov/fms-slab-ocean-model-technical-documentation/)
- [Bonan 2019 book](https://www.cambridge.org/us/academic/subjects/earth-and-environmental-science/climatology-and-climate-change/climate-change-and-terrestrial-ecosystem-modeling?format=HB&isbn=9781107043787)
- [Torres et al 19](https://link.springer.com/article/10.1007/s00382-018-4236-x?shared-article-renderer)
- [Thatcher & Jablonowski (2016)](https://gmd.copernicus.org/articles/9/1263/2016/gmd-9-1263-2016.pdf)
