# **Boundary flux and state exchanges**

# Calculation & accumulation
- fluxes are calculated using SurfaceFluxes.jl or a custom calculation, which can be called from:
    - the atmos model BC
        - can be specified as a coupled BC (e.g. [sea breeze example](https://github.com/CliMA/ClimaCoupler.jl/blob/as/agu-seabreeze/experiments/ClimaCore/sea_breeze/coupledbc.jl))
        - :) - fluxes applied at every substage of atmos timesteps (BCs updated at a higher resolution)
        - :( - aesthetically unappealing, since passing land state into atmos 
        - :) - interface similar to standalone runs 
    - the coupler
        - coupler extracts states and updates the BCs using a calllback 
        - :) - don't need to import land state into atmos BC
        - :( - cannot be done at each substage of the most resolved model 
        - :( - coupler needs to intervene at each timestep
        - :) - requires storage of calculated / accumulated fields within the coupler only, so it may alleviate workloads/memory for the node running the atmos model
- the coupler is responsible for regridding and transforming fluxes and states to the relevant modes, depending on what (complexity of) flux calculation is chosen
- the above should be consistent with the [CliMA Airspace plan](https://www.overleaf.com/project/6169b2b29040a9c1d73e2e38). 

*Currently, we are using the first approach as it facilitates prototyping (since all flux calculations are contained within the custom BC), but when optimising for performance we will compare both approaches.* 

## State passing for full SVAT scheme flux calculation
- should we pass last timestep or some average of the last coupling timestep? TBD

## Domain partition
- this is done via regridding onto an overlaping mesh - see [regrigging](regridding.md)

# Plan
- example for optimally estimating atmos state (average or interpolation) for land state calculation (not fluxes) *[not high priority]*