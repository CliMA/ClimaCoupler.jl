# **Two Coupled 2D LESs**
(WIP)

# Kelvin Helmholtz
- demonstrate conservation across rigid lid interface
- viscous stresses 

## Formulation
- builds on the sea breeze example
- compressible NS equations:


Density:
$$ \frac{\partial \rho}{\partial t} + \nabla \cdot ({\rho \vec{u}})= 0 $$

Momentum (flux form):
$$ \frac{\partial \rho \vec{u}}{\partial t} + \nabla \cdot ({\rho \vec{u} 
\otimes \vec{u} + pI})= \nabla \cdot (\rho \tau) - \rho g $$

Total Energy (excludes grav PE, as in Kang):
$$ \frac{\partial \rho e_{tik}}{\partial t} + \nabla \cdot ((\rho e_{tik} + p )\vec{u}) =   \nabla \cdot (\rho \tau \vec{u}) + \nabla \cdot \Big(\kappa_m\nabla T\Big) - \rho g \cdot \vec{u}, $$
where $e_{tik} = p (\rho(\frac{c_p}{c_v}-1))$ see appemdix for link with the full total energy.

<span style="color:grey">
Total Energy (includes grav):

$$ \frac{\partial \rho e_{tot}}{\partial t} + \nabla \cdot ((\rho e_{tot} + p )\vec{u}) = \nabla \cdot (\kappa \rho \nabla h_{tot}), $$
where $h_{tot}$ is the total specific enthalpy given by internal and potential energy contributions. 

Potential temperature (backup):
$$ \frac{\partial \rho \theta}{\partial t} + \nabla \cdot (\rho \theta \vec{u}) = \nabla \cdot (\kappa \rho \nabla \theta) $$
</span>


### Viscous forces
- stress tensor:
$$
\tau = \rho^{-1} \mu \Big(\nabla \vec(u) + \nabla(\vec{u})^T - \frac{2}{3} \mathcal{I} \nabla \cdot \vec{u}\Big)
$$
with $\mu$ the dynamic viscosity in Pa s. 



- Smagorinsky Closure:
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

## Appendix - $E_{tic} \rightarrow E_{tot}$
- total energy, $E_{tot} = E_{tpe} + E_{ke} $
- where total potential energy:
$$
E_{tpe} = GPE + IE = \int \rho g z dz + \int \rho c_v T dz = \int \rho c_p T dz
$$
where we used the ideal gas law and $c_p = R+c_v$, thus
$$
E_{tpe} = \int p \frac{c_p}{c_p-c_v} dz = \int p(1-\frac{c_v}{c_p})^{-1} dz 
$$
- thus $e_{tot} \neq e_{tik}$ where 
$$
e_{tik} = \frac{p}{\rho} (\frac{c_p}{c_v}-1)^{-1}
$$
with $\frac{c_p}{c_v} = c^2/ RT$



# Thermal Convection
- conservation and performance testing

# Refs:
- [kang et al 22](https://arxiv.org/pdf/2202.11890.pdf)
- [kang et al 21](https://arxiv.org/pdf/2101.09263.pdf)

