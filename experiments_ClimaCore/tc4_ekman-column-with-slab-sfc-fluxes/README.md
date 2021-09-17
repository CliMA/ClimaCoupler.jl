# Ekman Column + Slab with SurfaceFluxes.jl

This is a prototype for testing the surface flux model, using the ClimaAtmos Ekamn column and a slab, similar to Test Case 2 in the(see [Coupler Design Docs](https://www.overleaf.com/project/610c13492c7d0e8d459e72b8) for details). We're using the Nishiwaza 2018 to formulate the flux calculations below, as in `SurfaceFluxes.jl`. 


# Coupled boundary conditions

$$
\frac{\partial \rho \theta}{\partial t} = - \frac{\partial}{\partial z} \rho \overline{w'\theta'} 
$$
$$
\frac{\partial u }{\partial t} = - \frac{\partial}{\partial z} \overline{u'w'} 
$$


Where we parameterise:
$$
\overline{u'w'} = u*u* \approx - \tau / \rho = - <u> g_{am} = cd_m|<u>| <u>,
$$
$$
\overline{w'\theta'}= u*\theta* \approx - H / (\rho c_p) =g_{ac} (<\theta> - \theta_{sfc}) = cd_h  |<u>| (<\theta> - \theta_{sfc}).
$$

 Note that the units of $g$ are m/s; this is a slightly different $g$ definition compared with Bonan, Ch 7, where his $g$ is multiplied by a molar density. Figure 1.10 is helpful in clarifying the different representations of  turbulent fluxes across communities.   
    
Numerical iterations solve for 
$$
u* = \frac{\kappa}{\log(\Delta z/ z_{z0m}) - \Psi_m(\Delta z/ L) + \frac{z_{0m}}{\Delta z} \Psi_m (z_{0m}/L) + R_{z0m}[\Psi_m(z_{0m}/L) - 1]} <u> = \sqrt{cd_m} <u> = \sqrt{g_{am} <u> / \rho}
$$
$$
\theta* = \frac{\kappa / Pr}{\log(\Delta z/ z_{0h}) - \Psi_h(\Delta z/ L) + \frac{z_{0h}}{\Delta z} \Psi_h (z_{0h}/L) + R_{z0m}[\Psi_h(z_{0h}/L) - 1]} (<\theta> - \theta_{sfc})  = \sqrt{cd_h^2 / cd_m} (<\theta> - \theta_{sfc}) = \sqrt{g_{ac}^2/({g_{am} \rho <u>)}} (<\theta> - \theta_{sfc}) 
$$
with 
$$
<u> = \frac{1}{\Delta z} \int_0^{\Delta z} u(z) dz \,\,\,\,\,\,\,\,\,\,\,\,\, <\theta> = \frac{1}{\Delta z} \int_0^{\Delta z} \theta(z) dz
$$
and

$$
R_{z0m} = 1 - z_{0m} / \Delta z \,\,\,\,\,\,\,\,\,\,\,\,\, R_{z0h} = 1 - z_{0h} / \Delta z
$$
The Obukhov length is defined as:
$$
L  =  - \frac{u*^3 \overline{\theta}}{\kappa g \overline{w'\theta'|_s}} \,\,\,\,\,\,\,\,\,\,\,\,\, 
$$

where $\overline{\theta}$ is the basic potential temperature (?)


## Adapter for `SurfaceFluxes.jl`
- `dummy_surface_fluxes.jl` extracts the required drag/transfer coefficients from `SurfaceFluxes.jl`, i.e.:
$$cd_m (z)= \frac{u*^2}{  u(z)^2} \,\,\,\,\,\,\,\,\,\,\,\,\, cd_h(z) = \frac{u*\theta * }{ u(z) \Delta \theta (z)}$$
- the coefficients are then used to calculate the fluxes, similarly to the bulk formula (which assumes that the coefficients are constant) 

## Note 
- TODO: long-term physically meaningful tests once Div Operator is fixed