# Slab Ocean

The ocean can be coarsely parameterized as a 0 dimensional slab covering the whole globe, loosely representing an interactive well mixed layer with no dynamics. 

Following [Hansen et al. [67]](https://www.cesm.ucar.edu/models/atm-cam/docs/description/node48.html#hansen84), the evolution of the prognostic variable, ocean temperature, can be written as a heat flux equation (W m$^{-2}$):

$$
\rho_o h_o c_o \frac{\partial T_o}{\partial t} = (1-A)F_{ao} + AF_{io} + (1-A)F_{fr} + Q 
$$


where 
$\rho_o$ = ocean density (const)
$h_o$  = annual mean mixed layer depth (prescribed, or can be geographically variable - e.g. specified specified from [Levitus [105]](https://www.cesm.ucar.edu/models/atm-cam/docs/description/node48.html#levitus82))
$c_o$ = ocean specific heat capacity (const)

$F_{ao}$ = atmos <-> ocean flux = R$_{SW}$ + R$_{LW}$ + LH + SH + Q
- $Q$ = parameterizes the horizontal heat transport and deep water exchange (prescribed or based on the net surface energy flux from an atmosphere-only run with prescribed SSTs)
- R$_{SW}$, R$_{LW}$ = absorbed incoming solar short-wave and outgoing long-wave radiation energy fluxes at the surface 
- LH and SH are the latent and sensible heat fluxes from ocean to atoms (evaluated using $T_o$)

$F_{io}$ = sea ice <-> ocean flux  
$F_{fr}$ = heat flux gained from sea ice growth, nonzero only if $T_o<T_f$ (freezing point) and Q > 0

## Implementation
Using explicit time stepping the formulation is implemented as follows:
1) Top level adjustment (assuming it happens faster than deep adjustment)
	- $$T_o^{n+1} = T_o^n + 
	\frac{(1-A^n)F_{ao}^n}{\rho_o c_o h_o}\Delta t$$ 
2) adjust Q flux where $T_o^{n+1}  < 0$ K and $Q>0$ to avoid excessive cooling
	- $Q^n = Q^n\frac{T_f-T_o^{n+1}}{T_f}$ where $T_f = -1.8$C
	- additional $Q$ flux adjustments (e.g. tuning to observations) are possible in needed
3) Update $T_o^{n+1}$ with the new Q and ice-ocean fluxes
	- $T_o^{n+1} = T_o^{n+1} - \frac{Q_n + A^n F^n_{oi}}{\rho_o c_o h_o}\Delta t$
4) (if sea ice) Calculate $F_{frz}$ and adjust $T_o$ and $Q$ to points in warm-ocean waters ($T_o>0$C) 
	- $F_{frz}^{n+1} = \rho_o c_o h_o \max(T_f -T_o^{n+1} ,0) /\Delta t$
	- $T_o^{n+1} = \max(T_o^{n+1}, T_f)$
5) Normalise the Q flux to conserve energy
	- $Q = Q + (\overline{Q_0} - \overline{Q})(A_o/A_w)$ where bar = global average, $A_o$ = area of all ocean, $A_w$ = area of warm ocean(>0C) and Q is the unadjusted Q flux before step 2), so that $\overline{Q} = \overline{Q_0}$
	


## Values
- $\rho_o = 1.026e3$ km m-3
- $c_o = 3.93e3$ J km-1 K-1 
- $T_f$ = 1.8C
- 

## Rough parameterization of ocean surface fluxes 
- LW (upward pointing)
	- black body emission
	- LW flux = $F_{sfc}^u - F_{sfc}^d = \epsilon k T_o^4 - \epsilon*F_{atmos}$
		- $\epsilon = 0.98$ = emissivity
		- $k = 5.67e-8$ is the S. Boltzmann constant 
- SW (downward pointing)
	- SW flux = $F_{sfc}^d - F_{sfc}^u = \tau F_s - \alpha(\tau F_s )$ 
		- $F_s$ = solar flux
		- $\tau$ = transmissivity
		- $\alpha$ = albedo
- SHF = $c_p c_s\rho_a |V_0|(T_{o,skin} - T)$
	- $c_s = 1e-3 (1.2e-3)$ heat exchange coefficient (or over land) 
- LHR =  $c_s\rho_a |V_0|(q_{sat}(T_o,p_{sfc}) - q)$
	- $c_s= 1e-3$ heat exchange coeff over sea
References:
- https://www.cesm.ucar.edu/models/atm-cam/docs/description/node29.html



