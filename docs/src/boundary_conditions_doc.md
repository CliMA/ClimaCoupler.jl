# Boundary Conditions (BCs)

This section attempts to bridge mathematical considerations of BCs and their implementations in climate models.

## Types of boundary conditions
1. Dirichlet
    $$θ = f$$
2. Neumann
    $$\partial θ / \partial n = f$$
3. Robin
    $$c_1θ +  c_2\partial θ / \partial n = f$$
5. Cauchy
    $$θ = f_1$$  $$\partial θ / \partial n = f_2$$
4. Mixed
    - a mix of the above, each applied at a different section of the boundary

## Well-posedness and general context

- Well-posedness, or the availability of a unique robust solution, is easier to prove in some problems than others. For example, using the heat-diffusion equation this section demonstrates that some combinations of the above boundary conditions do not produce a unique solution. The problem of imposing boundary conditions in climate models is more complex, so well-posedness is more difficult to prove analytically, as discussed below. 

- Examples of imposing BCs:  

    1. **Diffusion equation in steady state**
        - The prognostic heat equation rewritten in terms of potential temperature, $θ$, is:
        $$\frac{\partial θ}{\partial t} = - ∇ \cdot (κ ∇ θ)$$
        - for simplicity, we assume a steady state ($\partial / \partial t = 0$), a diffusivity of $\kappa = 1$ m$^2$ s$^{-1}$ and that we're only in 1D (i.e. $∇ = \partial / \partial x$). Solving this requires two BCs in space (assuming we have no constraints in time), so we have: 
        $$\frac{\partial^2 θ}{\partial x^2} = 0$$
        - upon integration in $x$ we expect 
        $\frac{\partial θ}{\partial x} = c_1$ and $x = c_1x + c_2$, where $c_1$ and $c_2$ are set by the boundary conditions. 
        - If we have two Neumann boundary conditions (e.g. $dT/dx = 1$ at both $x_0$ = 1 and $x_N$ ), we can only obtain $c_1$ and the solution is not unique. We therefore need one of the boundaries to be constrained using a Dirichlet BC. 

    2. **Advection-diffusion equation with a linear flux**
        - stepping up the ladder of complexity is the advection diffusion problem:
            $$\frac{\partial θ}{\partial t} = - ∇ \cdot (\vec{u} \theta + κ ∇ θ) = - ∇ \cdot (F_{non-diffusive} + F_{diffusive})$$
        - [Miyaoka et al 17](https://www.scielo.br/pdf/tema/v18n2/2179-8451-tema-18-02-00253.pdf) have shown that by taking advantage of mass conservation, it is possible to derive a general BC for the advection-diffusion problem, as long as the flux is linear in $\theta$ (e.g. $\vec{u}(x,y,z)$). Evidently, this is still insufficient for climate modelling purposes.

    3. **Navier Stokes and climate models**
        - well-posedness of the NS equations is still unknown (after all it is one of the Clay Mathematics Institute's Millennium Problems)
        - in climate modelling it is often considered sufficient to use the general rule of thumb and use as many BCs as we have derivatives, without explicitly proving that the problem is analytically well posed. 
        - we implicitly assume mass conservation and the divergence theorem like in the [Miyaoka et al 17](https://www.scielo.br/pdf/tema/v18n2/2179-8451-tema-18-02-00253.pdf) study above
        - in practice, for global atmospheric models there are different constraints on each of the two boundaries:
            - lower boundary of the atmosphere
                - **Impenetrable**: reflective BCs are imposed by setting $w=0$ and the tangential wind remains unchanged. This is the `FreeSlip()` BC for momentum, `Impenetrable()` for mass and `Insulating()` for energy/moisture/tracers which may be used for atmosphere-only idealized studies.  
                - **Prescribed**: it is common to prescribe a constant (or variable if coupled to another model) surface temperature. This allows calculation of a surface sensible and latent heat fluxes, which are then added to any other sources of heat at the surface (e.g. radiation fluxes, etc). Once all surface fluxes are summed up, they are then passed to the atmospheric model as fluxes that enter the prognostic equations at the boundaries. E.g. see the `Prescribed()` BC, and `CoupledPrimary()`, `CoupledSecondary()` BCs for more complex setups.
            - upper boundary of the atmosphere
                - **Impenetrable**: reflective BCs are imposed by setting $w=0$ and the tangential wind remains unchanged. This is the `FreeSlip()` BC for momentum, `Impenetrable()` for mass and `Insulating()` for energy/moisture/tracers which may be used for atmosphere-only idealized studies.  
                - **Open**: it is undesirable to produce reflection in the upper levels, but open BCs are more difficult to implement
                - **Sponge layer**: an alternative to the above can be obtained by artificial (often linear or 2nd order) damping in the upper levels. The down side is that this sponge layer ofen takes up a substantial part of the model's depth.  

## Implementation in the ClimateMachine.jl

    Boundaries in DG:
        __________-  bc +
       | element  |  |  |            - = internal state
       | in       |  |  |           bc = boundary state 
       | question |  |  |            + = ghost state
       |__________|  |  |
                  \__ __/
                     V
        same location in space

- BCs in ClimateMachine.jl are imposed separately for the different compute kernels as illustrated by the [CMFlowChart](https://github.com/CliMA/ExperimentsMachine/blob/main/process_targetted_tests/hyperdiffusion/CM_flow.pdf).
- For the AtmosModel, the BCs that need to be considered are:


### Mass
- `Impenetrable()`: 
    - no normal component of mass flux, which means $u$ at $bc$:
        $u_{bc} = \vec{u}^- - (\hat{n} \cdot u^-) \hat{n}$
    - In conjunction with the *CentralNumericalFlux* at the boundary, this can be imposed via the ghost state as:

            # use CentralNumericalFlux with this reflective ghost state  ​
            function atmos_momentum_boundary_state!(_...)
                    ​state⁺.ρu = state⁻.ρu - 2 * dot(state⁻.ρu, n) .* SVector(n) 
            ​end
- `Penetrable()`: 
    - this means no boundary condition (free surface)

            # use the transmissive boundary flux: 
            function atmos_momentum_boundary_state!(_...)
                    ​state⁺.ρu = state⁻.ρu
            ​end


### Momentum 
- `NoSlip()`: 
    - this means $u$ vanishes are the boundary (i.e. Dirichlet BC), so that $u_{bc} = 0$
    - In conjunction with the *CentralNumericalFlux* at the boundary, this can be imposed via the ghost state as:

            function atmos_momentum_boundary_state!(nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient}, _...)
                state⁺.ρu = -state⁻.ρu
            end
            
            # a stabilising penalty term can be used here:
            numerical_boundary_flux_first_order!(_...) = ... end

            atmos_momentum_normal_boundary_flux_second_order!(_...) = nothing # no contribution from nF_diffusive

- `FreeSlip()`: $\nabla_h \cdot\vec{u} = 0$, $\hat{n} \cdot u = 0$
    - no diffusive normal flux at the boundary and no drag on the tangential components
    - $u_{bc}$ is determined by the mass flux BC above


            function atmos_momentum_boundary_state!(nf::NumericalFluxFirstOrder, _...)
                    state⁺.ρu -= 2 * dot(state⁻.ρu, n) .* SVector(n) # reflective to impose no normal flux
            end
            function atmos_momentum_boundary_state!(nf::NumericalFluxGradient, _...)
                    state⁺.ρu -= dot(state⁻.ρu, n) .* SVector(n) # non-reflective to capture the sign of the gradient
            end
            atmos_momentum_normal_boundary_flux_second_order!(_...) = nothing  # no contribution from nF_diffusive

- `DragLaw()`:
    - use momentum constraints on $u_{bc}$, just like FreeSlip(), but with a non-zero second-order normal flux (nF_diffusive) 

            function atmos_momentum_boundary_state!(nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient}, _...)
                atmos_momentum_boundary_state!(nf, Impenetrable(FreeSlip()), _...)
            end
            function atmos_momentum_normal_boundary_flux_second_order!(_...)
                fluxᵀn.ρu += C * state⁻.ρ * |u_h| * u_h # contribution from nF_diffusive
            end

### Energy / Moisture / Tracers
- `Insulating()`:
    - no normal diffusive flux (nF_diffusive) across the boundary (homogeneous Neumann BC) imposed via the second order flux, with the first-order and gradient fluxes assuming $u_{bc}$ from the mass and momentum BCs for consistency

            atmos_energy_boundary_state!(_...) 
                state⁺.ρe = state⁻.ρe
            end
            atmos_energy_normal_boundary_flux_second_order!( _...,) = nothing

- `Prescribed()`, `CoupledPrimary()`, `CoupledSecondary()`
    - same as Insulating(), but with a non-zero second-order normal flux (nF_diffusive)

            atmos_energy_boundary_state!(_...) 
                state⁺.ρe = state⁻.ρe
            end
            function atmos_energy_normal_boundary_flux_second_order!(_...)
                fluxᵀn.energy.ρe -= bc_energy.fn(state⁻, aux⁻, t)
            end

        where `bc_energy.fn` can be a flux value or a relaxation function to a surface temperature (e.g. see the slab land/ocean, bulk formulation, or the M-O flux-based NishizawaEnergyFlux function). NB: Dirichlet BC is not normally applied in this context.

## BC stability in DG for idealized problems (probably delete or expand?)
- boundary fluxes (numerical + physical) are implementable in DG via boundary numerical fluxes, between the interior point ($x^-$) and its exterior ghost equivalent ($x^+$). These need to be applied depending on the problem at hand. For example in the advective problem, the BCs will need to consider the flow direction, so that information is carried in the correct direction (i.e. downstream). 
- **Advection problem**
    - Imposed BCs
        - A) Impose BCs via numerical fluxes
            - e.g. central boundary flux (or another numerical flux, e.g. Rusanov):
            $$u\theta_{bc} = \frac{u\theta^- + u\theta^+}{2} = \frac{u\theta^- + f_0}{2}$$
        - B) Impose the flux directly at the boundary:
            $$u\theta_{bc} = f_0$$ 
            which is implemented via the ghost point as $\theta^+ = \theta^- + 2\theta_0 \rightarrow$  reflection principle
        - For an advective problem Lyapunov function analysis shows that A) is more dissipative (and stable) but B) is closer to the analytical conservation solution 
    - Free endpoints / open boundaries
        $$u\theta_{bc} = u\theta^-$$
        which is implemented as $\theta^+ = \theta^-  \rightarrow$ transmissive / free boundary (endpoint is just what it is)
- **Diffusion problem**
    - call the diffusive flux $\sigma = \kappa\partial_x \theta$
    - A) Neumann on $\theta$
        - is basically Dirichlet on $\sigma$, so that the boundary is transmissive $\theta$
        $$\sigma_{bc} = \frac{\sigma^- + \sigma^+}{2}$$
        and $\sigma^+ = \sigma^- + 2f_0$
        
        - $\theta_{bc} = \theta_0$ or $\sigma_{bc} = 0$ are stable and conservative
    - B) Dirichlet
        - transmissive on $\sigma_{bc}$ and explicit BCs on $\theta_{bc}$ are stable and conservative
- **Advection-diffusion problem**
    - best choice for the numerical solution to mimic the analytical solution is  $u\theta_{bc} = 0$ and $\sigma_{bc}=\sigma^-$
    - best choice for sufficient dissipation (and thus a more stable solution) is  $u\theta_{bc} = u\theta^- / 2$ and $\theta_{bc}=0$ and $\sigma_{bc}=\sigma^-$


## References:
- CliMA's numerics_old DesignDoc
- [useful Q&As](https://www.realclimate.org/index.php/archives/2009/01/faq-on-climate-models-part-ii/) 
- [Todd Lane's talk for climate modelling context](https://climateextremes.org.au/wp-content/uploads/2019/06/Introduction-to-Atmospheic-Modelling-1-Todd-Lane.pdf)