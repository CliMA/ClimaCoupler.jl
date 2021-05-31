# Slab Land

For prototyping of the coupler, we parameterize the land as a 0-dimensional slab covering the whole globe, analogous to the slab ocean with no water.

Analogous to the slab ocean of [Hansen et al. [67]](https://www.cesm.ucar.edu/models/atm-cam/docs/description/node48.html#hansen84), the evolution of the prognostic variable, surface soil temperature, can be written as a heat flux equation (W m$^{-2}$):

  
$$
\rho_s h_s c_s \frac{\partial T_{sfc}}{\partial t} = F_{tot} = R_{SW} - R_{LW} - SH + G - (LH)
$$

where
- R$_{SW} = F_{sfc}^d - F_{sfc}^u = \tau F_{sol} - \alpha(\tau F_{sol} )$ is the downward-pointing short-wave radiative flux
	- $F_{sol} = 1361$W m$^{–2}$ is the incoming solar radiation
	- $\tau = 0.9$ is the selected transmissivity of the atmosphere
	- $\alpha = 0.5$ is the selected albedo
	- $\rho_s =1500$ kg m$^{–3}$ is the selected density of soil (assumed constant)
	- $h_s= 2$ m is the selected soil depth
	- $c_s = 800$ J kg$^{–1}$ K$^{–1}$is the selected specific heat capacity of soil

- R$_{LW} = F_{sfc}^u - F_{sfc}^d = \epsilon \sigma \theta_{sfc}^4 - \epsilon F_{a}$ is the upward-pointing net long-wave flux
	- $\epsilon = 0.98$ is the selected emissivity or absorptivity of the surface
	- $\sigma = 5.67e^{-8}$ is the Stefan Boltzmann constant
	- $F_{a} =0$ = choosing no incoming longwave radiation from the atmosphere 

- SHF = $c_p g_a (T_{sfc} - T_a)$
	- $g_a$ = aerodynamic conductance (in complex models dependent on both land and atmosphere parameters, but in the single stack it is assumed constant, 2 mol m$^{-2}$ s$^{-1} \approx$ 58kg m$^{-2}$ s$^{-1}$ )
	- $c_p$ = specific heat capacity of the atmosphere at constant pressure (1004 J kg$^{–1}$ K$^{-1}$)

- $G = \kappa_s (\theta_{sfc}-\theta_h)/h_s$ is the soil storage
	- $\kappa_s = 2$ W m$^{-1}$ K$^{-1}$ is the selected thermal conductivity for soil (Fourier's law)
	- $\theta_h = 280$ K is the selected $\theta$ at $h_s$

- LHR = $\lambda g_w (q_{sat}(T_{sfc},p_{sfc}) - q_a)$
	- $g_w^{-1}=g_a^{-1}+g_s^{-1}$ = depends on aerodynamic conductance and surface conductance
	- $\lambda$ latent heat of vaporisation
	- $p_{sfc}$ is the surface pressure
	- $q_{sat}$ is the saturation specific humidity at the surface
	- $q_a$ is the specific humidity in the near surface atmospheric layer

The simple representation of the radiation will later be substituted by the radiation module in Atmos and a better representation of surface variables will be supplied by Atmos's SurfaceFluxes.jl (estimates based on the MO similarity theory) and subsequently by the Atmos's EDFM module. 

## Implementation in single stack

For this simple implementation, we assume a dry setup and that the aerodynamic conductance is constant. Using explicit time stepping the formulation is implemented as follows:

1) 

2) Top level adjustment (assuming it happens faster than deep adjustment)

$$T_s^{n+1} = T_s^n +\frac{F_{tot}^n}{\rho_s c_s h_s}\Delta t$$

3) pass $T_s^{n+1}$ back

## Implementation in LES (using MO)
1) initial Ts, x,  are passed to Atmos SurfaceFluxes.jl
    - using MO similary theory fluxes (e.g. $HF = u*\theta* * cp * \rho$) at the first Atmos nodal level are calculated and passed as the fTn condition
    - at each timestep within the coupling cycle, the total flux is calculated and accumulated in F_accum
    - during the flux calculation, atmos params chsnge while land params stay fixed (??)
2) At the end of atmos nsteps, coupling with land means sending F_accum to calculate T_sfc in the SlabLand
3) T_sfc is then passed to Atmos ans the cycle repeats

Other issues
- need to store F_accum externally
= e.g. when 3 models; need to reset after extracting, not before accumulation
- uis gc being updated at each timestep?
  

## References:
- Bonan book
- https://www.cesm.ucar.edu/models/atm-cam/docs/description/node29.html
- http://www.met.reading.ac.uk/~swrhgnrj/teaching/MT23E/mt23e_notes.pdf
