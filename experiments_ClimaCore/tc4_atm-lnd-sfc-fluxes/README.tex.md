# Atmos-Ocean-Soil Columns

This is a prototype for testing the surface flux model, using the atmos and land models of Test Case 3 in the [Coupler Design Docs](https://www.overleaf.com/project/610c13492c7d0e8d459e72b8) for details. We're using the Nishiwaza 2018 to formulate the flux calculations below, as in `SurfaceFluxes.jl`. 

# Coupled boundary conditions

$$
\frac{\partial \rho \theta}{\partial t} = - \frac{\partial}{\partial z} \rho \overline{w'\theta'} 
$$
$$
\frac{\partial u }{\partial t} = - \frac{\partial}{\partial z} \overline{u'w'} 
$$


Where we parameterise:
$$
\overline{u'w'} = u*u* \approx - \tau / \rho = - <u> g_{am} / \rho = cd_m|<u>| <u>,
$$
$$
\overline{w'\theta'}= u*\theta* \approx - H / (\rho c_p) =g_{ac} (<\theta> - \theta_{sfc}) / \rho = cd_h  |<u>| (<\theta> - \theta_{sfc})
$$

where numerical iterations solve for 
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
L  =  - \frac{u*^3 \overline{\theta}}{\kappa g \overline{w'\theta'|_s}} \,\,\,\,\,\,\,\,\,\,\,\,\, R_{z0m} = 1 - z_{0m} / \Delta z \,\,\,\,\,\,\,\,\,\,\,\,\, R_{z0h} = 1 - z_{0h} / \Delta z
$$
where $\overline{\theta}$ is the basic potential temperature (?)


## Variables supplied to SurfaceFluxes:
- in dynamic atmos mode (for FD):
    - $<u>$: $u$ averaged over first atmos layer
    - $\Delta z$: thickness of first atmos layer
- standalone land mode:
    - prescribe $u(Z)$ (instead of $<u>$)
    - $\Delta z = Z-d$, 

## Exchange variables 
- from land: $\theta_{sfc}$, $z_{0h}$, $z_{0m}$,
- from atmos: $<u>$, $<\theta>$, $<\rho>$, Pr 
- from `SurfaceFLuxes.jl`: 
$$cd_m = \frac{u*^2}{  u(z)^2} \,\,\,\,\,\,\,\,\,\,\,\,\, cd_h = \frac{u*\theta * }{ u(z) \Delta \theta (z)}$$
$$$$ 

 things depending on state, but can be assumed to be fixed over atmos step: T_land, relative_humidity land, for g_soil needs moisture in land

## Add to SurfaceFluxes.jl
- $g_{ac}$ (and $g_{s}$) calculation

## Other notes
- Land conductances added together for evaporation (this will be addresses later):
$$
g_{evap} = \frac{1}{ g_{ac}^{-1} + g_{sfc}^{-1}}
$$

