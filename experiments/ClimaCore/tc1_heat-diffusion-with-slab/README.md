# Test Case 1: Diffusive column with a slab

This is a basic test case, stripped of most interface wrappers for prototyping and tutorial purposes. 

## Model 1: ClimaCore diffusive column

The heat diffusion equation is 
$$
\partial_t T = \frac{\partial}{\partial z} \mu \frac{\partial}{\partial z}  T
$$
with boundary conditions of 
$\partial_t T =  - \frac{\partial}{\partial z}  F_{sfc}$   at the bottom of the column and
$T = 280 K$ at the top, and isothermal initial conditions. $\mu$ is a constant diffusion coefficient. 
We also use this model to calculate  fluxes, $F_{sfc}$, as described below.

## Model 2: Slab
The slab solves for temperature in a single layer, whose tendency is the accumulated fluxes divided by the coupling timestep plus a parameterisation of the internal processes, $G$ (in the default case we set this to zero).
$$
h_{lnd}  \, \partial_t T_{sfc} =  - F_{integ}  / \Delta t_{coupler}  + G
$$
where $h_{lnd}$ is a constant depth of the idealized slab layer. 

## Flux calculation and accumulation
This is done in Model 1 (which is assumed to have the shorter timestep). The calculation of (downward) surface fluxes, $F_{sfc}$, is:

$$
F_{sfc} = -\lambda (T_{sfc} - T)
$$

The fluxes are accumulated as a second 'prognostic' equation
$$
\partial_t F_{integ} =  - F_{sfc}.
$$

## Conservation tests
For this simple example, assume that the thermal heat capacity and density of both domains is one, so that the temperatures are analogous to specific energies in the conservation tests. Results can be viewed [here](https://www.overleaf.com/read/bgfmhgtncpws).

## TODO
- deploy CI + update conservation tests

