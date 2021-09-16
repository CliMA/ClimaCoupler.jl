# Columns with SurfaceFluxes.jl

This is a prototype for testing the surface flux model, using the atmos and land models of Test Case 3 in the [Coupler Design Docs](https://www.overleaf.com/project/610c13492c7d0e8d459e72b8) for details. We're using the Nishiwaza 2018 to formulate the flux calculations below, as in `SurfaceFluxes.jl`. 

There are two test cases:
- Heat diffusion minimal example:
    - this is an uncoupled integration test for surface fluxes (port to SurfaceFluxes when it's being revamped); test against bulk formula
- Coupling example: Atmos Ekman column + land column exchanging more realistic fluxes/states; flux calculation / accumulation is done by the coupler. 


# Coupled boundary conditions

$$
\frac{\partial \rho \theta}{\partial t} = - \frac{\partial}{\partial z} \rho \overline{w'\theta'} 
$$
$$
\frac{\partial u }{\partial t} = - \frac{\partial}{\partial z} \overline{u'w'} 
$$


Where we parameterise:
$$
\overline{u'w'} = u*u* \approx - \tau / \rho = - <u> g_{am} = cd_m|<u>| <u>,
$$
$$
\overline{w'\theta'}= u*\theta* \approx - H / (\rho c_p) =g_{ac} (<\theta> - \theta_{sfc}) = cd_h  |<u>| (<\theta> - \theta_{sfc}).
$$

 Note that the units of $g$ are m/s; this is a slightly different $g$ definition compared with Bonan, Ch 7, where his $g$ is multiplied by a molar density. Figure 1.10 is helpful in clarifying the different representations of  turbulent fluxes across communities.   
    
Numerical iterations solve for 
$$
u* = \frac{\kappa}{\log(\Delta z/ z_{z0m}) - \Psi_m(\Delta z/ L) + \frac{z_{0m}}{\Delta z} \Psi_m (z_{0m}/L) + R_{z0m}[\Psi_m(z_{0m}/L) - 1]} <u> = \sqrt{cd_m} <u> = \sqrt{g_{am} <u> / \rho}
$$
$$
\theta* = \frac{\kappa / Pr}{\log(\Delta z/ z_{0h}) - \Psi_h(\Delta z/ L) + \frac{z_{0h}}{\Delta z} \Psi_h (z_{0h}/L) + R_{z0m}[\Psi_h(z_{0h}/L) - 1]} (<\theta> - \theta_{sfc})  = \sqrt{cd_h^2 / cd_m} (<\theta> - \theta_{sfc}) = \sqrt{g_{ac}^2/({g_{am} \rho <u>)}} (<\theta> - \theta_{sfc}) 
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
L  =  - \frac{u*^3 \overline{\theta}}{\kappa g \overline{w'\theta'|_s}} \,\,\,\,\,\,\,\,\,\,\,\,\, 
$$

where $\overline{\theta}$ is the basic potential temperature (?)


## Variables supplied to SurfaceFluxes:
- in dynamic atmos mode (for FD; coupled to land or not):
    - $<u>$: $u$ averaged over first atmos layer
    - $\Delta z$: thickness of first atmos layer
- standalone land mode, driven by e.g. measurements of u at a given height $h$, $\theta$ at the same height $h$:
    - replace $<u>$ with $u(h)$
    - replace $\Delta z = h-d$, etc.

## Exchange variables 
- from land: $\theta_{sfc}$ (determined via prognostic land variables); $z_{0h}$, $z_{0m}$, $d$ (parameters of surface, at most slowly varying in time)
- from atmos: $<u>$, $<\theta>$, $<\rho>$, Pr 
- from `SurfaceFLuxes.jl`: 
$$cd_m (z)= \frac{u*^2}{  u(z)^2} \,\,\,\,\,\,\,\,\,\,\,\,\, cd_h(z) = \frac{u*\theta * }{ u(z) \Delta \theta (z)}$$
$$$$ 

 things depending on state, but can be assumed to be fixed over atmos step: T_land, relative_humidity land, for g_soil needs moisture in land

## Algorithm 
- this is called using the main `SurfaceFLuxes.jl` function
```
function surface_conditions(
    param_set::AbstractEarthParameterSet,
    MO_param_guess::AbstractVector,
    x_in::AbstractVector,
    x_s::AbstractVector,
    z_0::Union{AbstractVector, FT},
    θ_scale::FT,
    z_in::FT,
    scheme,
    wθ_flux_star::Union{Nothing, FT} = nothing,
    universal_func::Union{Nothing, F} = Businger,
    sol_type::NS.SolutionType = NS.CompactSolution(),
    tol::NS.AbstractTolerance = NS.ResidualTolerance{FT}(sqrt(eps(FT))),
    maxiter::Int = 10_000,
) where {FT <: AbstractFloat, AbstractEarthParameterSet, F}
```
- Define the function to calculate the `x - [monin_obukhov_length, u_star, theta_star]` that will be passed to the Newton solver
```
f!(F, x_all) = surface_fluxes_f!(F, x_all, args)
```
- Define Newton solver
```
nls = NS.NewtonsMethodAD(f!, MO_param_guess)
```
- Solve with the Newton solver from CliMA's `NonlinearSolvers.jl`
```
sol = NS.solve!(nls, sol_type, tol, maxiter)
```
using the general formula $x_{n+1} = x_n - \frac{f(x_n)}{\partial_x f(x_n)}$ 
```
function solve!(
    ::NewtonsMethodAD,
    x0::AT,
    x1::AT,
    f!::F!,
    F::FA,
    J::JA,
    J⁻¹::J⁻¹A,
    soltype::SolutionType,
    tol::AbstractTolerance{FT},
    maxiters::Int,
) where {FA, J⁻¹A, JA, F! <: Function, AT, FT}

    x_history = init_history(soltype, AT)
    F_history = init_history(soltype, AT)
    if soltype isa VerboseSolution
        f!(F, x0)
        ForwardDiff.jacobian!(J, f!, F, x0)
        push_history!(x_history, x0, soltype)
        push_history!(F_history, F, soltype)
    end
    for i in 1:maxiters
        f!(F, x0)
        ForwardDiff.jacobian!(J, f!, F, x0)
        x1 .= x0 .- J \ F
        push_history!(x_history, x1, soltype)
        push_history!(F_history, F, soltype)
        if tol(x0, x1, F)
            return SolutionResults(
                soltype,
                x1,
                true,
                F,
                i,
                x_history,
                F_history,
            )
        end
        x0 = x1
    end
    return SolutionResults(
        soltype,
        x0,
        false,
        F,
        maxiters,
        x_history,
        F_history,
    )
end
```

- 
## TODO minimal example (~TC1)
- note changes in SF and NS
- test: compare with bulk formula test
- note if still conservation problems due to divergence operator
- clean up @shows + document the SF functions (+ make README consistent)

## TODO physical example (TC3 minus ocean)
- add coupler flux calculation (also change this for TC3)
- couple stress, moisture
- add to SF: $g_{ac}$ (and $g_{s}$) calculation. $g_s$ needed only for evaporation and latent heat fluxes
- hook up to CI

## Other notes
- Land conductances added together for evaporation (this will be addresses later):
$$
g_{evap} = \frac{1}{ g_{ac}^{-1} + g_{s}^{-1}}
$$
- SHF also has a contribution due to evaporation (see Land Design Doc, Ch 2.8). So while the calculation of $u*$, etc, can be carried out as described, the total exchanged flux is not simply proportional to e.g. $u*\theta*$ for sensible heat. 

