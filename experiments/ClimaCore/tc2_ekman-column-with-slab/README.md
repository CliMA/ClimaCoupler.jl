# Test Case 2: Ekman column with a slab

The Ekamn column includes advection, diffusion and Coriolis terms. The slab land remains the same as in Test Case 1. This slightly more complex case extends the previous case, thereby enabling:
- using prognostic variables similar and equations to CliMA (can start thinking about BCs)
- using `ClimaAtmos.jl's` interface which will eventually be cast over all CliMA models
- zeroth-order inclusion of surface fluxes, which will eventually be replaced by importing `SurfaceFluxes.jl` 

Parameter descriptions and values for the equations below can be found [here](run.jl).

## Model 1: Ekman Column
The equations of motion for this problem are (see also section 5 in the Numerics chapter of the CliMADesignDocs):
$$
\partial_t \rho = - \partial_z (w \rho)
$$
$$
\partial_t \rho \theta= \partial_z ( \mu  \partial_z \theta - w \partial_z  \rho \theta)
$$
$$
\partial_t u =  \partial_z (\mu  \partial_z u) + f  (v - v_g) - w \partial_z u )
$$
$$
\partial_t v =  \partial_z (\mu  \partial_z v) - f  (u - u_g)  - w \partial_z v )
$$

The boundary conditions are impermeable for mass at both boundaries
$$ w\rho = 0.$$

For momentum and temperature, the top boundary is impenetrable and insulating (i.e. zero fluxes) and the bottom is forced by a diffusive flux:
$$\mu \partial_z u= F_u,$$   
$$\mu \partial_z v= F_v,$$   
$$\mu \partial_\rho\theta v= F_{\rho\theta},$$   
and at the top boundary the velocities are set equal to their geostrophic counterparts ($u = u_g$ and
$v = v_g$). All advective fluxes are set to zero at both boundaries. 

The vertical velocity tendency is set to zero at the boundaries where we solve: 
$$
\partial_t w_1 = - c_p \theta_1 \partial_z \Pi - \partial_z (\mu  \partial_z w) - g - w \partial_{z1} w
$$
where the subscript $_1$ refers to variables at the boundary faces, and the Exner function is $\Pi = (p_{\theta} / p_0)^{R / c_p}$.

We also use this model to calculate surface fluxes, $F_{sfc}$, as described below.

## Model 2: Slab
The slab solves for temperature in a single layer, whose tendency is the accumulated fluxes divided by the coupling timestep plus a parameterisation of the internal processes, $G$.
$$
\rho c h_s  \, \partial_t T_{sfc} =  - F_{integ}  / \Delta t_{coupler}  + G
$$
with
$$G = \kappa_s (T_{sfc}-T_h)/h_s$$
   where $T_h$ is the temperature at depth $h_s$ and the $_1$ subscript refers to variables extrapolated to the surface. For this simple example, we set the thermal capacity ($c$) and density ($\rho$) of both domains to be one. $G=0$ in the default example with conservation checks.


## Flux calculation and accumulation
We use Model 1 to calculate the thermal surface fluxes, $F_{sfc}$ (usually this should be done by the model with the shortest timestep). We parameterize the surface fluxes as:
$$
F_{sfc} = R_{SW} - R_{LW} - SH 
$$
The zeroth order representation of the short (SW) and long (LW) wave radiation fluxes and sensible heat flux (SH) are calculated as:  
    $$SH = cp_d g_a (\theta_{sfc} - \theta_a) = cp_d g_a (\theta_{sfc} - \theta_1)$$
    $$R_{SW} = F_{sfc}^d - F_{sfc}^u = (1-\mu)\tau F_{sol} (1 + \sin(2 \pi \tau_d^{-1} t))$$
    $$R_{LW} = F_{sfc}^u - F_{sfc}^d = \epsilon \sigma \theta_{sfc}^4 - \epsilon F_{a}$$
      
The fluxes are accumulated as anogther second 'prognostic' equation of Model 1 as
$$
\partial_t F_{integ} =  - F_{sfc}.
$$

## Conservation tests
For this simple example, assume that the thermal heat capacity and density of both domains is one, so that the temperatures are analogous to specific energies in the conservation tests. Results can be viewed [here](https://www.overleaf.com/read/bgfmhgtncpws).
