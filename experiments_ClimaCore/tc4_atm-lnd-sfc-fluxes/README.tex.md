# Atmos-Ocean-Soil Columns

This is a prototype for testing the three-way model coupling. Please Test Case 4 in the [Coupler Design Docs](https://www.overleaf.com/project/610c13492c7d0e8d459e72b8) for details. 

# Coupled boundary conditions

$$
\frac{\partial \rho \theta}{\partial t} = - \frac{\partial}{\partial z} \rho \overline{w'\theta'} 
$$
$$
\frac{\partial u }{\partial t} = - \frac{\partial}{\partial z} \overline{u'w'} 
$$


Where we parametrerise:
$$
\overline{u'w'} = u*u*
$$
$$
\overline{w'\theta'} = u*\theta*
$$

where numerical iterations solve for 
$$
u* = \frac{\kappa}{\log(\Delta z/ z_{z0m}) - \Psi_m(\Delta z/ L) + \frac{z_{0m}}{\Delta z} \Psi_m (z_{0m}/L) + R_{z0m}[\Psi_m(z_{0m}/L) - 1]} <u> = \sqrt{cd_m} <u>
$$
with 
$$
<u> = \frac{1}{\Delta z} \int_0^{\Delta z} u(z) dz
$$

similar for $\theta*$.

Land: conductances added together for evaporation

$$
\overline{w'\theta'} = c_H |u| (\theta_{sfc} - T_{atmos} ) 
$$

$(u \times g) = C_h\times u^2= \rho * u*^2$ # for momentum
$H = cp*(theta - T_land) *g = rho c_p u* theta* = c_p*(theta - T_land)*c_h*u?$
"c_h_eff * u" = 1/((C_h*u)^-1 + g_soil^-1)?

in dynamic atmos mode = <u>, \Delta z
in driving land mode , know u(Z) = uZ, set \Delta z = Z-d, <u> = uZ


Exchange of parameters between models
- e.g. roughness_length_map(x,y,z,t) - 
- 
 things depending on state, but can be assumed to be fixed over atmos step: T_land, relative_humidity land, for g_soil needs moisture in land