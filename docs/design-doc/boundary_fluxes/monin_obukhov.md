# **Monin Obukhov Theory Application for Boundary conditions**

# General formulation

We need to parameterise the SGS portion of turbulent sensible, latent and heat fluxes:
$$
\frac{\partial \rho \theta}{\partial t} = ... - \frac{\partial}{\partial z} \rho \overline{w'\theta'} 
$$
$$
\frac{\partial u }{\partial t} = ... - \frac{\partial}{\partial z} \overline{u'w'} 
$$
$$
\frac{\partial q}{\partial t} = ... - \frac{\partial}{\partial z}  \overline{w'q'} 
$$

where the primes denote perturbations from the average environmental values. The scales of these turbulent motions are:
$$
\overline{u'w'} \approx  -u*u* \,\,\,\,\,\,\,\,\,\, \overline{w'\theta'} \approx  -u*\theta* \,\,\,\,\,\,\,\,\,\, \overline{w'q'} \approx  -u*q*
$$


Where the turbulent fluxes can be written in various formulations:

$$
\rho \overline{u'w'} = - \tau  \approx \rho K_m \frac{\partial \overline{u}}{\partial z} \approx - g_{am}  (<u> - u_{sfc})  = \rho C_D |<u>|  (<u> - u_{sfc}),
$$
$$
\rho \overline{w'\theta'}= - SH / (c_p) \approx \rho K_h \frac{\partial \overline{\theta}}{\partial z} \approx g_{ac} (<\theta> - \theta_{sfc}) = \rho C_H  |<u>| (<\theta> - \theta_{sfc})
$$
$$
\rho \overline{w'q'}= - LH / (L) \approx \rho K_h \frac{\partial \overline{q}}{\partial z} \approx g_{ac} (<q> - q*_{sfc}) = \rho C_H  |<u>| (<q> - q*_{sfc})
$$
This theory therefore advocates that we could estimate turbulent eddy fluxes from mean state variables and an estimation of their respective efficiencies, $C_D$ and $C_H$. These efficiencies can be assumed constant (reculting in a *bulk aerodynamic formula*) or they can be estimated using more sophisticated methods. One of such methods is the Monin Obukhov theory, which (with later improvements for representing near-surface variables in FV) assumes:
    $$
    u* = \frac{\kappa}{\log(\Delta z/ z_{z0m}) - \Psi_m(\Delta z/ L_{MO}) + \frac{z_{0m}}{\Delta z} \Psi_m (z_{0m}/L_{MO}) + R_{z0m}[\Psi_m(z_{0m}/L_{MO}) - 1]} <u> 
    $$
    $$
    \theta* = \frac{\kappa / Pr}{\log(\Delta z/ z_{0h}) - \Psi_h(\Delta z/ L_{MO}) + \frac{z_{0h}}{\Delta z} \Psi_h (z_{0h}/L_{MO}) + R_{z0m}[\Psi_h(z_{0h}/L_{MO}) - 1]} (<\theta> - \theta_{sfc})  
    $$
    $$
    q* = \frac{\kappa / Pr}{\log(\Delta z/ z_{0h}) - \Psi_h(\Delta z/ L_{MO}) + \frac{z_{0h}}{\Delta z} \Psi_h (z_{0h}/L_{MO}) + R_{z0m}[\Psi_h(z_{0h}/L_{MO}) - 1]} (<q> - q*_{sfc})
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
    L_{MO}  =  - \frac{u*^3 \overline{\theta}}{\kappa g \overline{w'\theta'|_s}} =  - \frac{u*^3}{\kappa \overline{w'b'|_s}} = \frac{u*^2}{\kappa b*} \,\,\,\,\,\,\,\,\,\,\,\,\, 
    $$
    where $\overline{\theta}$ is the average potential temperature and where:
    $$
    b^* =\frac{g}{\rho} \frac{\kappa / Pr}{\log(\Delta z/ z_{0h}) - \Psi_h(\Delta z/ L_{MO}) + \frac{z_{0h}}{\Delta z} \Psi_h (z_{0h}/L_{MO}) + R_{z0m}[\Psi_h(z_{0h}/L_{MO}) - 1]}((<\rho> - \rho_{sfc})
    $$
    The FV-corrected universal functions $\Psi_m$ and $\Psi_h$ which are derived semi-empirically (see NK19).


# Implementation
- The implementation in CliMA is in SurfaceFluxes.jl and is split into two stages:
- 1. solving $L_{MO}$
    - If we substitute the velocity and buoyancy scales into the definition of $L_{MO}$, given $z_1$ (first interior point) and $z_{sfc}$, the only unknown parameter is $L_{MO}$. No easy solution exists, so this equation needs to be solved numerically until a certain tolerance level is reached. We use the Newton iterator.
    - input are roughness lengths and variables at surface and first interior model level:
        
        
            #Roughness lengths
            z0m = FT(0.001)
            z0b = FT(0.001)
            state_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
            state_in = SF.InteriorValues(z_in, u_in, ts_in)

            kwargs = (; state_in, state_sfc, z0m, z0b)
            sc = SF.ValuesOnly{FT}(; kwargs...)
            uf = UF.Businger() # default
            result = SF.surface_conditions(param_set, sc, uf)


- 2. conversion to efficiencies:
    $$
    C_D  = \frac{u*^2}{(<u> - u_{sfc})^2}  \,\,\,\,\,\,\,\,\,\,\,\,\, C_H  = \frac{u*\theta*}{(<\theta> - \theta_{sfc})^2}  \,\,\,\,\,\,\,\,\,\,\,\,\,
    $$
    - to avoid singularity in $C_D$ a gustiness parameter may be employed. In stable conditions (singularity in $C_H$, neutral conditions) there is no buoyancy flux and the momentum coefficient follows the log law of a wall and is a function of z:
    $$
    C_D = \Big(\frac{\kappa}{\ln{\Delta z / z_{ob}}}\Big)^2
    $$ 
- 3. knowing $L_{MO}$ and scales allows us to compute the scales, and in turn the profile of a variable between $z_{sfc}<z<z_1$, somply by rearranging the definitions for $u*$, $\theta*$ and $q*$ for $u(z)$, $\theta(z)$ and $q(z)$

            function recover_profile(
                param_set::APS,
                sc::AbstractSurfaceConditions,
                L_MO,
                Z,
                X_in,
                X_sfc,
                transport,
                uft::UF.AUFT,
                scheme::Union{FVScheme, FDScheme},
            ) 

- 4. collect output
    - Result struct of type SurfaceFluxConditions{FT} contains:
        - L_MO:   Monin-Obukhov lengthscale
        - shf:    Sensible Heat Flux
        - lhf:    Latent Heat Flux
        - ρτxz:   Momentum Flux (Eastward component)
        - ρτyz:   Momentum Flux (Northward component)
        - ustar:  Friction velocity
        - Cd:     Momentum Exchange Coefficient
        - Ch:     Thermal Exchange Coefficient

# Refs
- Nishizawa & Kitamura 18 (NK18)
- [SurfaceFluxesDesignDoc](https://www.overleaf.com/project/616d983d2ac263356dc204dd)
- [CanopyAirspaceDesignDoc](https://www.overleaf.com/project/6169b2b29040a9c1d73e2e38)