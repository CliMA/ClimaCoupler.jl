#=
# Overview

This sea breeze simulation consists of an atmosphere above ocean and land slabs. The difference in heating between the land and ocean components drives circulation: cool ocean air flows towards the land at the surface while warm air over land rises and flows over the ocean. The physical models are described below.

## Flux calculation
$$F_{sfc} = c_p \rho_1 C_H |u_1| (\theta_{sfc} - \theta_{atm1})$$
$$F_{integ} = \int_{\Delta t_{coupler}} F_{sfc}  dt$$
where $c_p$ is the specific heat, $C_H = 0.0015$ is the bulk transfer coefficient for sensible heat, $|u_1|$ is the near-surface atmospheric wind speed. $F_{integ}$ has units of $J m^{-2}$. $\theta_{sfc}$ represents the potential temperature at the land or ocean surface.

We assume that the potential temperature is defined with respect to the surface pressure, so that $\theta_{sfc} = T_{sfc}$. 

## Problem Domain 

## Initial, Forcing and Boundary Conditions

- Atmospheric boundary conditions:
    - Horizontal: periodic
    - External vertical: free slip + insulating
    - Coupled vertical: $F_{sfc}$ boundary flux

- Initial conditions

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

# References
- [Antonelli & Rotunno 2007](https://journals.ametsoc.org/view/journals/atsc/64/12/2007jas2261.1.xml?tab_body=pdf)
=#
