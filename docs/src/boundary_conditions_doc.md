# Boundary Conditions (BCs)

This section is an attempt to bridge mathematical considerations of BCs and their implementations in climate models.

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
    - a mix of the above, each applied at a different part of the boundary

## Well-posedness and general context
- well-posedness, or the availability of a unique robust solution, is easier to prove in some problems than others. For example, for the heat diffusion equation, the next section shows that some combinations of the above boundary conditions do not produce a unique solution. The problem in climate models is, however, more complex, and we discuss it in the section after next. 

    1. BCs for diffusion equation in steady state
        - The prognostic heat equation rewritten for temperature $θ$ is:
        $$\frac{\partial θ}{\partial t} = - ∇ \cdot (κ ∇ θ)$$
        - for simplicity, we assume a steady state ($\partial / \partial t = 0$), $\kappa = 1$ m$^2$ s$^{-1}$ and that we're only in 1D (i.e. $∇ = \partial / \partial x$). Solving this required 2 BCs in space (assuming we have no constraints in time), so we have: 
        $$\frac{\partial^2 θ}{\partial x^2} = 0$$
        - upon integration in $x$ we expect 
        $\frac{\partial θ}{\partial x} = c_1$ and $x = c_1x + c_2$, where $c_1$ and $c_2$ are set by the boundary conditions. 
        - If we have 2 Neumann boundary conditions (e.g. $dT/dx = 1$ at both $x_0$ = 1 and $x_N$ ), we can only obtain $c_1$ and the solution is not unique. We therefore need one of the boundaries to have a Dirichlet BC. 


    2. BCs for advection-diffusion equation with a linear flux
        - stepping up the complexity is the advection diffusion problem:
            $$\frac{\partial θ}{\partial t} = - ∇ \cdot (\vec{u} \theta + κ ∇ θ) = - ∇ \cdot (F_{non-diffusive} + F_{diffusive})$$
        - [Miyaoka et al 17](https://www.scielo.br/pdf/tema/v18n2/2179-8451-tema-18-02-00253.pdf) have shown that by taking advantage of mass conservation, it is possible to derive a general BC for the advection-diffusion problem, as long as the flux is linear in $\theta$ (e.g. $\vec{u}(x,y,z)$). Evidently, this is still insufficient for climate modelling purposes.

    3. Navier Stokes and Climate Models
        - well posedness of the NS equations is still unknown (after all it is one of the Clay Mathematics Institute's Millennium Problems)
        - in climate modelling it is often considered sufficient to use the general rule of thumb and use as many BCs as we have derivatives, without explicitly showing that the problem is well posed. 
        - we implicitly assume mass conservation and the divergence theorem like in the study above
        - in practice, for GCMs there are different constraints on the two boundaries
            - lower boundary of the atmosphere
                - **Impenetrable**: reflective BCs are imposed by setting $w=0$ and the tangential wind unchanged. This is the `FreeSlip()` BC for momentum, `Impenetrable()` for mass and `Insulating()` for energy/moisture/tracers which may be used for atmosphere-only idealized studies.  
                - **Prescribed**: it is common to prescribe a constant (or variable if coupled to another model) surface temperature. This allows calculation of a surface sensible heat flux, which is then added to any other sources of heat at the surface (e.g. radiation, latent heat flux, etc). Once all fluxes are calculated they are then passed to the model as basically a **Robin boundary condition**, and is a sum of all diffusive and nondiffusive fluxes). E.g. see the `Prescribed()` BC, and `CoupledPrimary()`, `CoupledSecondary()` BCs for more complex setups.
            - upper boundary of the atmosphere
                - **Impenetrable**: reflective BCs are imposed by setting $w=0$ and the tangential wind unchanged. This is the `FreeSlip()` BC for momentum, `Impenetrable()` for mass and `Insulating()` for energy/moisture/tracers which may be used for atmosphere-only idealized studies.  
                - **Open BCs**: it is undesirable to produce reflection in the upper levels, but open BCs are more difficult to implement
                - **Sponge layer**: an alternative to the above can be obtained by artificial (often linear or 2nd order) damping in the upper levels 

## Implementation in the ClimateMachine.jl

    Boundaries in DG:
                 -  bc +
        element  |  |  |            - = internal state
        in       |  |  |           bc = boundary state 
        question |  |  |            + = ghost state
                 \__ __/
                    V
        same location in space
### Mass
- `Impenetrable()`: 
    - no normal component of mass flux, which means $u$ at $bc$:
        $u_{bc} = \vec{u}^- - (\hat{n} \cdot u^-) \hat{n}$
    - In conjunction with *CentralNumericalFlux* at the boundary, this can be imposed via the ghost state as:

            # use CentralNumericalFlux with this reflective ghost state  ​
            function atmos_momentum_boundary_state!(_...)
                    ​state⁺.ρu = state⁻.ρu - 2 * dot(state⁻.ρu, n) .* SVector(n) 
            ​end
- `Penetrable()`: 
    - this means no bounary condition (free surface)

            # use the transmissive boundary flux: 
            function atmos_momentum_boundary_state!(_...)
                    ​state⁺.ρu = state⁻.ρu
            ​end


### Momentum 
- `NoSlip()`: 
    - this means u vanishes are the boundary (i.e. Dirichlet BC), so that $u_{bc} = 0$
    - Inconjunction with *CentralNumericalFlux* at the boundary, this can be imposed via the ghost state as:

            function atmos_momentum_boundary_state!(
                nf::NumericalFluxFirstOrder, _...)
                state⁺.ρu = -state⁻.ρu
            end
            function atmos_momentum_boundary_state!(
                nf::NumericalFluxGradient, _...)
                state⁺.ρu = zero(state⁺.ρu) # DesignDocs say state⁺.ρu = -state⁻.ρu
            end
            
            # this can be used for a stabilising penalty term:
            numerical_boundary_flux_first_order!(_...) = ... end

            atmos_momentum_normal_boundary_flux_second_order!(_...) = nothing

- `FreeSlip()`: $\nabla_h \cdot\vec{u} = 0$, $\hat{n} \cdot u = 0$

        function atmos_momentum_boundary_state!(
            nf::NumericalFluxFirstOrder
                state⁺.ρu -= 2 * dot(state⁻.ρu, n) .* SVector(n)
        end
        function atmos_momentum_boundary_state!(
            nf::NumericalFluxGradient,
                state⁺.ρu -= dot(state⁻.ρu, n) .* SVector(n)
        end
        atmos_momentum_normal_boundary_flux_second_order!(_...) = nothing

- `DragLaw()`:

        function atmos_momentum_boundary_state!(
            nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient}, _...)
            atmos_momentum_boundary_state!(nf, Impenetrable(FreeSlip()), _...)
        end
        function atmos_momentum_normal_boundary_flux_second_order!(_...)
            fluxᵀn.ρu += C * state⁻.ρ * |u_h| * u_h
        end

### Energy / Moisture / Tracers
- `Insulating()`:

        atmos_energy_boundary_state!(_...) = nothing
        atmos_energy_normal_boundary_flux_second_order!( _...,) = nothing
- `Prescribed()`, `CoupledPrimary()`, `CoupledSecondary()`

        atmos_energy_boundary_state!(_...,) = nothing

        function atmos_energy_normal_boundary_flux_second_order!(_...)
            fluxᵀn.energy.ρe -= bc_energy.fn(state⁻, aux⁻, t)
        end
    where `bc_energy.fn` can be a flux value or a relaxation function to a surface temperature (e.g. see slab land/ocean, or bulk formulation, or the M-O flux-based NishizawaEnergyFlux function). NB: Dirichlet is not normally applied in this context.

## BC stability in DG for idealized problems (probably delete)
- boundary fluxes (numerical + physical) are implementable in DG via boundary numerical fluxes, between the interior point ($x^-$) and its exterior ghost equivalent ($x^+$). These need to be applied depending on the problem at hand. If using an advective problem, the BCs will need to consider the direction, so that the flow carries information in the right direction. 
- Advection fluxes
    - Imposed BCs
        - A) Impose BCs via numerical fluxes
            - e.g. central boundary flux (or another numerical flux, e.g. Rusanov):
            $$u\theta = \frac{u\theta^- + u\theta^+}{2} = \frac{u\theta^- + f_0}{2}$$
        - B) Impose directly the flux at the boundary:
            $$u\theta = f_0$$ 
            so that $\theta^+ = \theta^- + 2\theta_0$ = reflection principle
        - For an advective problem Lyapunov function analysis shows that A) is more dissipative (and stable) but B) is closer to the analytical conservation solution 
    - Free endpoints / open boundaries
        $$u\theta = u\theta^-$$
        so that $\theta^+ = \theta^-$ = transmissive / free boundary (endpoint is just what it is)
- Diffusive fluxes
    - call the diffusive flux $\sigma = \kappa\partial_x \theta$
    - A) Neumann on $\theta$
        - is basically Dirichlet on $\sigma$, so that the boundary is transmissive $\theta$
        $$\sigma = \frac{\sigma^- + \sigma^+}{2}$$
        and $\sigma^+ = \sigma^- + 2f_0$
        
        - $\theta = \theta^-$ or $\sigma = 0$ are stable and conservative
    - B) Dirichlet
        - transmissive on $\sigma$ and explicit BCs on $\theta$
- Advective and diffusive
    - best choice for the numerical solution to mimic the algebraic solution is  $u\theta = 0$ and $\sigma=\sigma^-$
    - best choice for dissipation (and more stable solution) is  $u\theta = u\theta^1 / 2$ and $\theta=0$ and $\sigma=\sigma^-$


## References:
- [Bonan 2019 book](https://www.cambridge.org/us/academic/subjects/earth-and-environmental-science/climatology-and-climate-change/climate-change-and-terrestrial-ecosystem-modeling?format=HB&isbn=9781107043787)
- https://www.cesm.ucar.edu/models/atm-cam/docs/description/node29.html
- http://www.met.reading.ac.uk/~swrhgnrj/teaching/MT23E/mt23e_notes.pdf
- useful Q&As: https://www.realclimate.org/index.php/archives/2009/01/faq-on-climate-models-part-ii/ 
- https://climateextremes.org.au/wp-content/uploads/2019/06/Introduction-to-Atmospheic-Modelling-1-Todd-Lane.pdf
- old clima numerics design docs