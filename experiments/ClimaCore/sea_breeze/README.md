


## Atmosphere Conservation Equations

Density:
$$ \frac{\partial \rho}{\partial t} + \nabla \cdot ({\rho \vec{u}})= S(\chi, ...) $$

Momentum (flux form):
$$ \frac{\partial \rho \vec{u}}{\partial t} + \nabla \cdot ({\rho \vec{u} \otimes \vec{u} + pI})= \nabla \cdot (\rho \tau) - \rho g + F_{B}(...)$$


Potential temperature:
$$ \frac{\partial \rho \theta}{\partial t} + \nabla \cdot (\rho \theta \vec{u}) = \nabla \cdot (\kappa \rho \nabla \theta) $$


Total Energy (possibly replace potential temperature equation with total energy conservation):
$$ \frac{\partial \rho e_{tot}}{\partial t} + \nabla \cdot ((\rho e_{tot} + p )\vec{u}) = \nabla \cdot (\kappa \rho \nabla h_{tot}), $$
 
where $h_{tot}$ is the total specific enthalpy given by internal and potential energy contributions. 

Tracer transport: 
$$ \frac{\partial \rho \chi}{\partial t} + \nabla \cdot (\rho \chi \vec{u}) = \nabla \cdot (\kappa \rho \nabla \chi) + S(\chi, ...)$$

Smagorinsky Closure:
$$ 
\rho\tau = -2\rho\nu\vec{S} 
$$
$$ 
\vec{S} = \frac{1}{2}((\nabla u) + (\nabla u)^{T})
$$
and 
$$ 
\nu = (C_{s}\Delta_{x,y,z})^2\sqrt{2S_{ij}S_{ij}}
$$

with $\Delta_{x,y,z}$ the grid lengthscale (sometimes approximated as a geometric average $\Delta = (\Delta_x\Delta_y\Delta_z)^{1/3}$), $\nu$ is the kinematic viscosity (calculated here with the Smagorinsky model), $\vec{S}$ the symmetric rate-of-strain tensor, $\tau$ the diffusive momentum flux tensor. 

## Slab land ODE
$$\rho_l c_l H_l \partial_t T_{lnd} = - F_{integ} / \Delta t_{coupler}$$
- where $\rho_l = 1500$ kg m $^{-3}$, $c_l=800$ J K $^{-1}$ kg $^{-1}$, $H_l=1$ m are the density, specific heat and depth of the land slab. 
- $F_{integ}$ are the integrated surface fluxes in time (see below)

## Slab ocean ODE
$$\rho_o c_o H_o \partial_t T_{ocn} = - F_{integ} / \Delta t_{coupler}$$
- where $\rho_o = 1025$ kg m $^{-3}$, $c_o=3850$ J K $^{-1}$ kg $^{-1}$, $H_o = 100$ m are the density, specific heat and depth of the land ocean. 

## Flux calculation
$$ F_{sfc} = c_p \rho_1 C_H |u_1| (\theta_{sfc} - \theta_{atm1}) $$
$$ F_{integ} = \int_{\Delta t_{coupler}} F_{sfc}  dt $$
where $c_p$ is the specific heat, $C_H = 0.0015$ is the bulk transfer coefficient for sensible heat, $|u_1|$ is the near-surface atmospheric wind speed. $F_{integ}$ has units of $J m^{-2}$. $\theta_{sfc}$ represents the potential temperature at the land or ocean surface.

We assume that the potential temperature is defined with respect to the surface pressure, so that $\theta_{sfc} = T_{sfc}$. 
## Problem Domain 

3D, 
$$
L_{x}, L_{y}, L_{z} = (50km , 5km, 4km)
$$
$$
\Delta_{x}, \Delta_{y}, \Delta_{z} = (100m , 100m, 100m)
$$


## Initial, Forcing and Boundary Conditions

- Atmospheric boundary conditions:
    - Horizontal: periodic
    - External vertical: free slip + insulating
    - Coupled vertical: $F_{sfc}$ boundary flux

- Initial conditions

<!-- 
Stably stratified initial conditions with given buoyancy frequency $$N = \sqrt{2} \times 10^{-2} s^{-1}$$ for the control simulation.

Heat flux:
$$HF = \frac{1}{2}Q_L^*(\tanh((x - x_{0})/\lambda)+1)$$ 
Drag coefficient:
$$C_{D} = \frac{1}{2}C_D^*(\tanh(x - x_{0}/\lambda)+1)$$ 
$$\lambda = 100 ~\mathrm{m}$$
Geostrophic Wind forcing (momentum terms): 
$$
\vec{GW} = 
\rho f( \vec{u}- \vec{u_{g}})
$$

via equations (8), (9) in Antonelli and Rotunno. Control simulation with $$C_D^* = 0.007,$$  
$$Q_L^* = 0.078 \mathrm{~mKs^{-1}}.$$ 

$x_{0}$ is some offset coordinate based on the domain construction.  -->


##  Smagorinsky

Consider components of the viscous stress tensor in three dimensions
$$
\tau_{xx} = 2\nu \frac{\partial u}{\partial x},
$$
$$
\tau_{yy} = 2\nu \frac{\partial v}{\partial y},
$$
$$
\tau_{zz} = 2\nu \frac{\partial w}{\partial z} ,
$$
$$
\tau_{xy} = \nu \Big(\frac{\partial u}{\partial y} +  \frac{\partial v}{\partial x}\Big),
$$
$$
\tau_{xz} = \nu \Big(\frac{\partial u}{\partial z} +  \frac{\partial w}{\partial x}\Big),
$$
$$
\tau_{yz} = \nu \Big(\frac{\partial v}{\partial z} +  \frac{\partial w}{\partial y}\Big).
$$
Assume terms in the $y$-direction are neglected (2-dimensional simplicfication). The contributions to the momentum equation are then given by: 
$$
(\rho u):  \partial_{x} (\rho \tau_{xx}) + \partial_{z}(\rho\tau_{xz})  = \partial_x  \Big(2\nu \frac{\partial u}{\partial x}\Big) + \partial_z\Big(\nu \frac{\partial u}{\partial z}\Big) + \partial_z\Big(\nu \frac{\partial w}{\partial x}\Big),
$$
$$
(\rho w): \partial_{x} (\rho \tau_{zx})+ \partial_{z}(\rho\tau_{zz})  = \partial_x\Big(\nu \frac{\partial u}{\partial z}\Big) +  \partial_x\Big(\nu \frac{\partial w}{\partial x}\Big) + \partial_z\Big(2\nu\frac{\partial w}{\partial z} \Big). \\
$$
Which can be interpreted as, for horizontal-momentum:
1) Horizontal divergence of vertical gradients of cell-centered variables $u$
2) Vertical divergence of vertical gradients of cell-centered variables $u$
3) Vertical divergence of horizontal gradients of cell-face variables $w$

and for vertical-momentum, as:
1) Horizontal divergence of vertical gradients of cell-centered variables $u$ $TODO: Check Geometry.transform construction ???$
2) Horizontal divergence of horizontal gradients of cell-face variables $w$
3) Vertical divergence of vertical gradients of cell-face variables $w$

Thermal diffusivities are related to the modelled eddy viscosity through the turbulent Prandtl number which takes a typical value of $Pr_{t}= 1/3$ such that $\kappa = \nu/Pr_{t}$

## Required Operators 

Spectral (horizontal) divergence
$$
\mathrm{hdiv} = \mathrm{O.Divergence(~)}
$$


Spectral (horizontal) gradient 
$$
\mathrm{hgrad} = \mathrm{O.Gradient(~)}
$$


Finite Difference (vertical) gradient 
$$
\mathrm{hgrad} = \mathrm{O.GradientC2F(~)}
$$



# References
- [Antonelli & Rotunno 2007](https://journals.ametsoc.org/view/journals/atsc/64/12/2007jas2261.1.xml?tab_body=pdf)
