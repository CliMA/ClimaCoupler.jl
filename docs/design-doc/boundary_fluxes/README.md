# **Boundary flux and state exchanges**
- the coupler is responsible for regridding and transforming fluxes and states to the relevant models

# Turbulent fluxes

## Flux calculation
### Internal energy
- we use three approaches of varying complexity
    - simple
        $$
        F_{sfc} = - λ (T_{sfc} - T_v1) 
        $$
    - aerodynamic formulae + radiative flux
        - bulk
            $$
            F_{sfc} = R_{SW} - R_{LW} - SH - LH
            $$
            - where
                $$
                SH = g_a c_{pd} (\theta_{sfc} - \theta_1) = \rho_1 C_H c_{pd} |V|_1 (\theta_{sfc} - \theta_1) 
                $$
                $$
                LH = g_w L (q_{sat}(θ_{sfc}, p_{sfc}) - q_a) = \rho_1 C_L L (q_{sat}(θ_{sfc}, p_{sfc}) - q_a)
                $$
                - with the heat and moisture transfer coefficients $C_H = C_L \approx 0.0011$ , $|V| = \sqrt{(u^2+v^2)^2}$, specific dry-air heat capacity $c_{pd} = 1004 J K^{-1} kg^{-1}$, latent heat of vaporization $L=2.5e6 J K^{-1}$
                - saturation specific humidity, $q* = p_v*(T) / (\rho R T)$ (this may require saturation adjustment)
            
        - Monin Obukhov
            - same as above but calculate $C_H$ and $C_L$ (and $C_D$) using Monin Obukhov similarity theory
            - based on input roughness lengths, atmos near-surface state variables and surface state variables it calculates the turbulent fluxes and exchange coefficients. 
            - see [these notes](monin_obukhov.md) for more implementation details


### Moisture
- moisture follows similar approaches to internal energy           
    $$
    P = g_w (q_{sat}(θ_{sfc}, p_{sfc}) - q_a) = \rho_1 C_L (q_{sat}(θ_{sfc}, p_{sfc}) - q_a)
    $$

### Momentum
- aerodynamic formulae
    - bulk
        $$
        τ = \rho C_D * |V|_1 (u_1 - u_sfc) 
        $$
        - where $C_D\approx 0.0015$ is the drag coefficient 
    - Monin Obukhov
        - $C_D$ calculate as part of the method described above

## Flux accumulation
- flux accumulation is used the time integration of the finest model (e.g. atmos), which requires defining an additional "prognostic" variable: 
$$
d(F_{integrated})/dt  = F_{sfc}
$$
$$
F_{accumulated} = {F_{integrated}} / Δt_{coupler}
$$

# State passing / time interpolating
- decide when to pass last timestep vs computing some  running average of the last coupling timestep


# Implementation 

## Location of calculations
- fluxes are calculated using SurfaceFluxes.jl or a custom calculation, which can be called from:
    - the atmos model BC (implemented)
        - can be specified as a coupled BC (e.g. [sea breeze example](https://github.com/CliMA/ClimaCoupler.jl/blob/as/agu-seabreeze/experiments/ClimaCore/sea_breeze/coupledbc.jl))
        - :) - fluxes applied at every substage of atmos timesteps (BCs updated at a higher resolution)
        - :( - aesthetically unappealing, since passing land state into atmos 
        - :) - interface similar to standalone runs 
    - the coupler (to be implemented)
        - coupler extracts states and updates the BCs using a calllback 
        - :) - don't need to import land state into atmos BC
        - :( - cannot be done at each substage of the most resolved model 
        - :( - coupler needs to intervene at each timestep
        - :) - requires storage of calculated / accumulated fields within the coupler only, so it may alleviate workloads/memory for the node running the atmos model
        - :) - possibility to calculate fluxes on the exchange (finest) grid

    - *Currently, we are using the first approach as it facilitates prototyping (since all flux calculations are contained within the custom BC), but when optimising for performance we will compare both approaches.* 


## Domain partition
- using the linear remapping we can remap to split differently-resolved sub-domains, as long as their combined area = area of the whole domain
- this can be done via the linear mapping operator (and masking if necessary)
- see [regrigging](regridding/README.md)

# Supporting test cases
- dry (advection-)diffusion column + slab
    - see [TC1]("../../experiments/ClimaCore/tc1_heat-diffusion-with-slab/run.jl"), [TC2]("../../experiments/ClimaCore/tc2_ekman-column-with-slab/README.md")
    - simple coupling
    - field exchange
        - col > slab: 
            - heat flux
            - (momentum flux)
        - slab > col: 
            - $\theta_{sfc}$
- miost sea-breeze+slab
    - SurfaceFluxes.jl
    - uses ClimaSim BC
    - field exchange
        - col > slab: 
            - heat flux:
            - momentum flux:
            - moisture flux: 
        - slab > col: 
            - $\theta_{sfc}$
            - roughness
- convection 2D LES ([Kang et al 22](https://arxiv.org/pdf/2202.11890.pdf))
    - SurfaceFluxes.jl
    - uses ClimaSim BC
    - comparison with literature - performance
    - field exchange
        - LES1 > LES2: 
            - heat flux:
            - momentum flux:
            - moisture flux: 
        - LES2 > LES1: 
            - $\theta_{sfc}$
            - $q_{sfc}$
            - $u,v$
            - roughness
- RCE / Ekman column + land
    - accumulate states continuously
    - all equations for AMIP
    ...
- AMIP (FMS)
    - land < - > atm 

        - t_surf
        - u_surf
        - v_surf
        - q_surf

        - frac_open_sea
        - ex_ice_frac
        - land_frac
        - ex_rough_moist
        - ex_rough_heat
        - albedo (diff bands?)

        - ex_flux_sw
        - ex_flux_lw

        - ex_p_surf
        - wind
        - t
        - q
        - t_flux
        - q_flux
        - u_flux
        - v_flux
        - dtaudu
        - dtaudv

        - u_star
        - b_star
        - q_star
        - ex_wind
        - ex_cd_q
        - ex_cd_t
        - ex_cd_m

    - ocn /ice > atm
        - t_surf
        - u_surf
        - v_surf
        - q_surf


- CMIP

# Misc notes
- in FMS u and v stress components are calculated and regridded separately
- q* and z_0 also need to come from land  /ocean
- radiative heating of the air within the canopy airspace is ignored

# Questions
- shouldn't we use virtual T in C_H above?
- constant canopy height? 
- wind speed or wind vector necessary by canopy? 

# CliMA Strategy

## AMIP
- Ensure SurfaceFluxes complete
- 

## Optimisation
- example for optimally estimating atmos state (average or interpolation) for land state calculation (not fluxes) *[not high priority]*

# Refs
- [CliMA Airspace plan](https://www.overleaf.com/project/6169b2b29040a9c1d73e2e38). 