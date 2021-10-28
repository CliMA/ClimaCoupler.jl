


### Atmosphere Conservation Equations

Density:
$$ \frac{\partial \rho}{\partial t} + \nabla \cdot ({\rho \vec{u}})= S(\chi, ...) $$

Momentum:
$$ \frac{\partial \rho \vec{u}}{\partial t} + \nabla \cdot ({\rho \vec{u} \otimes \vec{u} + pI})= \nabla \cdot (\rho \tau) - \rho g + F_{B}(...)$$

Potential temperature:
$$ \frac{\partial \rho \theta}{\partial t} + \nabla \cdot (\rho \theta \vec{u}) = \nabla \cdot (\kappa \rho \nabla \theta) $$

Smagorinsky Closure:
$$ \frac{\partial \rho \chi}{\partial t} + \nabla \cdot (\rho \chi \vec{u}) = \nabla \cdot (\kappa \rho \nabla \chi) + S(\chi, ...)$$

$$ 
\rho\tau = -2\rho\nu\vec{S} 
$$
$$ 
\vec{S} = \frac{1}{2}((\nabla u) + (\nabla u)^{T})
$$
and 
$$ 
\nu = (C_{s}\Delta_{x,y,z})^2\sqrt{2S:S}
$$

with $\Delta_{x,y,z}$ the grid lengthscale (sometimes approximated as a geometric average), $\nu$ is the kinematic viscosity (calculated here with the Smagorinsky model), $\vec{S}$ the symmetric rate-of-strain tensor, $\tau$ the diffusive momentum flux tensor. 

### Slab land ODE
$$\partial_t T_{lnd} = - F_{integ} / \Delta t_{coupler}$$
- where $F_{integ}$ are the integrated surface fluxes in time:
$$ F_{integ} = \int -\lambda_{lnd} (T_{sfc} - T_{atm}) dt_{coupler} $$
where $\lambda_{lnd}$ absorbs the land properties relative to the atmos properties, and the drag coefficient: $c_D(\rho_a c_p)(\rho_l c_l h_l)^{-1}$

### Slab ocean ODE
$$\partial_t T_{ocn} = - F_{integ} / \Delta t_{coupler}$$
- where $F_{integ}$ are the integrated surface fluxes in time:
$$ F_{integ} = -\lambda_{ocn} (SST - T_{atm}) $$

### Initial, Forcing and Boundary Conditions

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

$x_{0}$ is some offset coordinate based on the domain construction. 



### Problem Domain 

3D, 
$$
L_{x}, L_{y}, L_{z} = (50km , 5km, 4km)
$$
$$
\Delta_{x}, \Delta_{y}, \Delta_{z} = (100m , 100m, 100m)
$$

### Required Operators 

Spectral (horizontal) divergence
$$
\mathrm{hdiv} = \mathrm{O.Divergence(~)}
$$


Spectral (horizontal) gradient 
$$
\mathrm{hgrad} = \mathrm{O.Gradient(~}
$$


Spectral (horizontal) gradient 
$$
\mathrm{hgrad} = \mathrm{O.Gradient(~}
$$