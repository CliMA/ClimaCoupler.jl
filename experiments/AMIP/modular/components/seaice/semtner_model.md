# **Simple Sea Ice Model**

# Semtner's Zero-layer Model

Sea ice is approximated as a model that does not occupy a particular gridpoint between the mixed layer (ocean) and the atmosphere. It is essentially implemented in each domain column as an ODE, with no horizontal transport between the columns. In the absence of ice, the sea-ice model reduces to the slab ocean formulation. The ice is assumed to have a negligible heat capacity (so there is no energy storage due to internal temperature changes of the ice). The only storage changes arise from the ice thickness changes, which result from temperature differences at the ice surface or at the ice base. 

We followed the FMS implementation as in [Zhang et al 22](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2021MS002671), which can be found on [GitHub](https://github.com/sally-xiyue/fms-idealized/blob/sea_ice_v1.0/exp/sea_ice/srcmods/mixed_layer.f90) and is itself a modification of the 0-layer model of [Semtner 1976](https://www.atmosp.physics.utoronto.ca/people/guido/PHY2502/articles/seaice-landice/AJSemtner_1976.pdf) (appendix). Chronologically, the algorithm follows these steps, with all fluxes defined positive upward:

## 1. Ice thickness, $h_i$
$$
L_i \frac{dh_i}{dt} = F_{atm} - F_{base}
$$
- with the latent heat of fusion, $L_i=3 \times 10^8$ J m$^{-3}$, and where the (upward-pointing) flux into the atmosphere is
$$
F_{atm} = F_{rad} + F_{SH} + F_{LH} \approx \lambda (-{T_{sfc}} - T_{atm})
$$
- where the latter approximation was used for testing in earlier prototypes. $F_{atm}$ will be obtained from the atmospheric model (via the coupler). The flux at the ice base from the mixed layer is
$$
F_{base} = F_0(T_{ml} - T_{melt})
$$
- where $T_{melt} = 273.16$ K is the freezing temperature, and the basal heat coefficient $F_0 = 120$ W m$^{-2}$ K$^{-1}$. 

## 2. Ocean mixed layer temperature, $T_{ml}$
- $T_{ml}$ is the standard slab ocean formulation in ice-free conditions:
$$
\rho_w c_w h_{ml}\frac{dT_{ml}}{dt} = - F_{atm}
$$
- while ice-covered conditions require that:
$$
\rho_w c_w h_{ml}\frac{dT_{ml}}{dt} = - F_{base}
$$

## 3. Transitions between ice free and ice covered conditions
- If the updated $T_{ml}^{t+1} < T_{melt}$, set $T_{ml}^{t+1} = T_{melt}$ and grow ice ($h_i^{t+1}$) due to the corresponding energy deficit.
- If the updated $h_i^{t+1} <= 0$ from a non-zero $h_i^t$, adjust $h_i^{t+1} = 0$ and use the surplus energy to warm the mixed layer. 

## 4. Surface temperature ($T_s$)
- $T_s$ is determined implicitly using a balance between $F_{atm}(T_s)$ and the conductive heat flux through the ice slab, $F_{ice}$:
$$
F_{atm} = F_{ice} = k_i \frac{T_{melt} - T_s}{h_i}
$$
- where $k_i = 2$ W m$^{-2}$ k$^{-1}$ is the thermal conductivity of ice.
- currently the implicit solve is implemented as one Newton iteration:
$$
T_s^{t+1} = T_s + \frac{F}{dF /d T_s}  = T_s^{t} + \frac{- F_{atm}^t + F_{ice}^{t+1}}{k_i/h_i^{t+1} + d F_{atm}^t / d T_s^t}   
$$
- where $h_i^{t+1}$ is the updated $h^i$ from the previous section, and $d F_{atm}^t / d T_s^t$ needs to be supplied from the atmosphere model (or crudely calculated in the coupler, given $T_s$, turbulent diffusivities and transfer coefficients, and atmos state). 
- Where $T_s^{t+1} > T_{melt}$, we set $T_s^{t+1} = T_{melt}$. Where there is no ice $T_s^{t+1} = T_{ml}^{t+1}$.

## 5. Update ice mask
- `mask = 1` if ice and `mask = 0` if no ice

# Q flux (optional)
- We can add an additional flux to the RHS of the $T_{ml}$ equations in (2), which corresponds to a more realistic ocean heating, coarsely mimicking otherwise neglected ocean dynamics, such as lateral advection, convection and diffusion. This is especially needed to improve low latitude oceanic forcing.

- An analytic formulation can be written as:
$$
Q = Q_0(1-2\phi^2/w_\phi^2) \frac{exp(- (\phi^2/w_\phi^2))}{cos(\phi)}
$$
- where $\phi$ is latitude in radians, $Q_0$ is the amplitude of the equatorial heating and $w_\phi$ the width of the heating in radians. 

# Alternatives

## [Semtner 1976](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2021MS002671) 3-layer model 
- Semtner (1976) developed a simple model for the evolution ice (2 layers) and snow (1 layer) temperatures, which is represented by a 1D diffusve process also accounting for radiation, melting and energy release from brine poskets and accumulating snow. 
- Main formation of sea ice - energy fluxes from vertical boundaries and heat storage in brine pockets. The vertical processes retain the central role.  
- Semtner 1976 found that the 0-layer model is broadly comparable in terms of accuracy as the 3-layer model

## Ocean-based dynamical sea-ice model
- to be coordinated after AMIP


