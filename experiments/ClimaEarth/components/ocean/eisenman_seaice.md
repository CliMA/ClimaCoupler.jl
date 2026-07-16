# Eisenman-Zhang Model
Thermodynamic 0-layer model, based on the Semtner 1976 model and later refined by
Eisenman & Wettlaufer (2009) and Zhang et al. (2021).

There are three prognostic variables, the height of ice (`h_i`), the ocean mixed layer depth (`T_ml`) and the surface air temperature (`T_s`).

Note that Eisenman sea ice assumes gray radiation, and no snow coverage.

## Formulation
In ice-covered conditions:
$$
L_i \frac{dh_i}{dt} = F_{atm} - F_{base} - Q
$$
with the density-weighted latent heat of fusion, $L_i=3 \times 10^8$ J m$^{-3}$, and where the (upward-pointing) flux into the atmosphere is
$$
F_{atm} = F_{rad} + F_{SH} + F_{LH}
$$
where $F\_{rad}$ is the net radiative flux, $F_{SH}$ is the sensible heat flux and $F_{LH}$ is the latent heat flux.

The basal flux is
$$
F_{base} = F_0(T_{ml} - T_{melt})
$$
where $T_{melt} = 273.16$ K is the freezing temperature, and the basal heat coefficient $F_0 = 120$ W m$^{-2}$ K$^{-1}$.

The lateral oceanic heat flux is parameterized as a prescribed $Q$-flux. With $F_{base} = T_{ml}^{t+1} = T_{melt}$ when ice is present, and zero $Q$-flux in the default setup, the base flux is zero in that case.

The surface temperature of the ice-covered surface is solved by balancing $F_{atm}(T_s)$ and the conductive heat flux through the ice slab, $F_{ice}$:
$$
F_{atm} = F_{ice} = k_i \frac{T_{melt} - T_s}{h_i}
$$
where $k_i = 2$ W m$^{-2}$ K$^{-1}$ is the thermal conductivity of ice.
Currently the solve is implemented as one Newton iteration (sufficient for the current spatial and temporal resolution - see Semtner, 1976):
$$
T_s^{t+1} = T_s + \frac{F}{dF /d T_s}  = T_s^{t} + \frac{- F_{atm}^t + F_{ice}^{t+1}}{k_i/h_i^{t+1} + d F_{atm}^t / d T_s^t}
$$
where $h_i^{t+1}$ is the updated $h^i$ from the previous section, and $d F_{atm}^t / d T_s^t$ needs to be estimated using $T_s + \delta T$.

### Warm surface
In ice-free conditions, the ocean temperature $T_{ml}$ assumes the standard slab model representation:
$$
\rho_w c_w h_{ml}\frac{dT_{ml}}{dt} = - F_{atm}
$$
In this case, $T_s^{t+1} = T_{ml}^{t+1}$

### Frazil ice formation
The frazil ice formation rate is parameterized as a function of the mixed layer temperature. It occurs when the newly calculated ocean temperature $T_{ml}^{t+1} < T_{melt}$. Since the mixed layer is not allowed to cool below $T_{melt}$, the energy deficit is used to grow ice:
$$
\frac{dh_i}{dt} = \Delta h_i^{t+1} + \frac{\Delta T \rho_w c_w h_{ml}}{L_i \Delta t}
$$
with $\Delta T =  T_{melt} - T_{ml}^{t+1} $.

### Transition to ice-free conditions
- If the updated $h_i^{t+1} < 0$ from a non-zero $h_i^t$, the ice height is set to zero and the surplus energy warms the mixed layer.

## Potential extensions
- add area `ice_area_fraction` adjustment (e.g., assuming a minimal thickness of sea ice, below which the grid area becomes part ice and part ocean)

# References
- [Semtner 1976](https://journals.ametsoc.org/view/journals/phoc/6/3/1520-0485_1976_006_0379_amfttg_2_0_co_2.xml)
- [Eisenman & Wettlaufer 2009](https://www.pnas.org/doi/full/10.1073/pnas.0806887106)
- [Zhang et al 2021](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2021MS002671), whose implementation can be found on [GitHub](https://github.com/sally-xiyue/fms-idealized/blob/sea_ice_v1.0/exp/sea_ice/srcmods/mixed_layer.f90)
