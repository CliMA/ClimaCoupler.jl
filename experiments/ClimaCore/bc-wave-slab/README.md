# **Baroclinic Wave + Slab**

# Atmosphere
The momentum equations are in the advective form, and tracers in the consevative form, namely:

- Density:
$$ \frac{\partial \rho}{\partial t} + \nabla \cdot ({\rho \vec{u}})= 0 $$

- Momentum (flux form):
$$ \frac{\partial \vec{u_h}}{\partial t} + \vec{u} \cdot \nabla \vec{u_h} = - \frac{1}{\rho}\nabla_h p
+ \frac{\partial}{\partial z} K_v \frac{\partial}{\partial z} \vec{u_h}
$$
$$ \frac{\partial w}{\partial t} + \vec{u} \cdot \nabla w= 
- \frac{1}{\rho}\frac{\partial p}{\partial z}
- \nabla_z \Phi 
+ \frac{\partial}{\partial z} K_v \frac{\partial}{\partial z} w 
$$

- Total energy:
$$ \frac{\partial \rho e_{tot}}{\partial t} + \nabla \cdot (\rho h_{tot} \vec{u}) = \frac{\partial}{\partial z} K_v \frac{\partial}{\partial z} h_{tot}
$$

where the total specific enthalpy and total specific  energy are
$$ h_{tot} =  e_{tot} + \frac{p}{\rho}  \,\,\,\,\,\,\,\, \,\,\,\,\,\,\,\, e_{tot} = c_v T + \Phi + \frac{1}{2}\vec{u}^2 
$$
(note that $h_{tot} \neq h = c_vT + p/\rho = c_p T$, the specific enthalpy in the thermodynamic sense), $\Phi = gz$ is the geopotential,
$u_h$ is the horizontal velocity vector, $w$ the vertical velocity, $\rho$ the density, $p$ pressure, $K_v$ the vertical diffusivity (assumed constant here). 

## Boundary conditions (BCs)
- We implement BCs similarly to other climate models. 
    - First-order fluxes (i.e., advective fluxes) are always set to zero, corresponding to the *free-slip* and *impenetrable* BC, where:
    $$
    w = 0 \,\,\,\,\,\,\, \partial_t w = 0 \,\,\,\,\,\,\, \nabla \times\vec{u_h}=0 \,\,\,\,\,\,\, \nabla \cdot \vec{\rho u_h}=0  \,\,\,\,\,\,\, \nabla \cdot \rho h_{tot} \vec{u_h}=0
    $$
    - Second-order fluxes (i.e., diffusive fluxes)
        - `NoFlux()`: By default we have *impenetrable* or *insulating* BCs (no second-order fluxes) at all boundaries. 
        - `BulkFormula()`: Applied to tracers (e.g., temperature and moisture), this imposes a boundary fluxes (e.g., sensible and latent heat) calculated using the bulk aerodynamic formulae  using prescribed surface values of ($T_{sfc}$ and $q_{sfc}^{sat}$).  At the surface, the bulk sensible heat flux formula for total enthalpy essentially replaces the above: 
            $$ (K_v \rho \partial_z h_{tot})_{sfc}$$
            For **total energy**, we have two choices:
            - 1. enthalpy flux: 
                $$ (K_v \rho \partial_z h_{tot})_{sfc} \rightarrow 
                \hat{n} \cdot  \rho C_H ||u||^{1} (h^1- h_{sfc})   
                = F_S
                $$
            - 2. sensible (and latent) heat flux. The sensible heat flux is: 
                $$ (K_v \rho \partial_z h_{tot})_{sfc} \rightarrow 
                \hat{n} \cdot C_H c_{pd} ρ^{1} ||u||^{1} (T^{1} - T_{sfc}) 
                +  \hat{n} \cdot C_H ρ^{1} ||u||^{1} (\Phi^{1} - \Phi_{sfc})      
                = F_S
                $$
                where $^{1}$ corresponds to the lowest model level, $C_H$ is the dimensionless thermal transfer coefficient, $c_{pd}$ is the specific heat capacity for dry air $||u||$ the wind speed. This is the *bulk turbulent sensible heat flux* parameterization, and $F_S$ is positive when atmosphere receives energy from the surface. 
                The contribution of the kinetic energy is usually O(1e4) smaller and is neglected, but it can be added to F_S as:
                $$
                F_{S_{tot}} = F_S + \hat{n} \cdot C_D ρ^{1} ||u||^{1} (\vec{u_h}^{1})^2 
                $$        
         - ` DragLaw()`: essentially the bulk formula for momentum       
            $$ \frac{\partial}{\partial z} K_v \frac{\partial}{\partial z} \vec{u_h} \rightarrow 
            \hat{n} \cdot C_D ρ^{1} ||u||^{1} \vec{u_h}^{1}          
            = F_M
            $$


        - `BulkFormulaCoupled()`: same as `BulkFormula()`, but $T_{sfc}$ is passed from the neighboring model.

    - The diffusive fluxes are applied via the `vertical_diffusion` ClimaAtmos model sub-component. To apply boundary fluxes without diffusion in the atmospheric interior, the viscosity coefficient needs to be set to zero: $ν = FT(0)$
- Current setup
    - total energy and momentum:
        - at z=0, we use `BulkFormulaCoupled()` with the enthalpy flux formulation, and the `DragLaw()`, with $C_D = C_H = 0.001$
        - interior diffusivity is set to $\nu = 5$ m^2/s
        - the the top $F_S = F_M = 0$
    - All other boundary fluxes are set to 0.

## Initial conditions
- we initialize with a perturbation in a balanced background state, as in:
https://climate.ucdavis.edu/pubs/UMJS2013QJRMS.pdf

# Heat Slab
The slab solves for temperature in a single layer, whose tendency is the accumulated fluxes divided by the coupling timestep plus a parameterisation of the internal processes, $G$.
$$
\rho c h_s  \, \partial_t T_{sfc} =  - F_{integ}  / \Delta t_{coupler}  
$$

# Conservation checks
- this uses the `sum` of `ClimaCore/Fields/mapreduce.jl`, which produces a sum weighted by the area Jacobian. 

# NB:
- first coupled iteration does not call rhs!s
- slab `T_sfc` gets huge numbers when using `SSPRK33`. ok with `Euler`

# TODO
- calculate fluxes at cell faces (not centers) - interp
- need better CA interface for specifying sensible / latent heat fluxes - i.e. fluxes in terms of temp and q, not enthaplpy fluxes
- add more precise flux accumulation
- conservation tests - add error threshold and exception, interval, show option, and make a general interface for it

# References
- [Kang et al 2021](https://arxiv.org/abs/2101.09263)




