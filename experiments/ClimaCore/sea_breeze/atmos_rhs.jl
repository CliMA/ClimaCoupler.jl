# # Atmospheric Model

#=
## Atmosphere Conservation Equations

Density:
```math
\frac{\partial \rho}{\partial t} + \nabla \cdot ({\rho \vec{u}})= S(\chi, ...).
```

Momentum (flux form):
```math
\frac{\partial \rho \vec{u}}{\partial t} + \nabla \cdot ({\rho \vec{u} \otimes \vec{u} + pI})= \nabla \cdot (\rho \tau) - \rho g + F_{B}(...).
```

Potential temperature:
```math
\frac{\partial \rho \theta}{\partial t} + \nabla \cdot (\rho \theta \vec{u}) = \nabla \cdot (\kappa \rho \nabla \theta).
```

Total Energy (possibly replace potential temperature equation with total energy conservation):
```math
\frac{\partial \rho e_{tot}}{\partial t} + \nabla \cdot ((\rho e_{tot} + p )\vec{u}) = \nabla \cdot (\kappa \rho \nabla h_{tot}),
```

where ``h_{tot}`` is the total specific enthalpy given by internal and potential energy contributions.

Tracer transport:
```math
\frac{\partial \rho \chi}{\partial t} + \nabla \cdot (\rho \chi \vec{u}) = \nabla \cdot (\kappa \rho \nabla \chi) + S(\chi, ...).
```

Diffusion (Constant Viscosity):
The simplest model to represent diffusive processes is a constant-viscosity model, with
prescribed kinematic viscosity ``\nu`` such that the stress tensor can be modelled by
```math
\rho\tau = -2\rho\nu\nabla u.
```

Smagorinsky Closure:
The Smagorinsky closure is an eddy-viscosity model that captures the effect of energy
transfer to the smallest scales of motion in the flow.
```math
\begin{aligned}
\rho\tau &= -2\rho\nu\vec{S}, \\
\vec{S} &= \frac{1}{2}((\nabla u) + (\nabla u)^{T}), \\
\nu &= (C_{s}\Delta_{x,y,z})^2\sqrt{2S_{ij}S_{ij}}.
\end{aligned}
```

with $\Delta_{x,y,z}$ the grid lengthscale (sometimes approximated as a geometric average
``\Delta = (\Delta_x\Delta_y\Delta_z)^{1/3}``), $\nu$ is a spatially varying kinematic viscosity
that depends on the local shear, ``\vec{S}`` the symmetric rate-of-strain tensor,
``\tau`` the diffusive momentum flux tensor. In stratified flows, we can apply a correction
to the eddy viscosity to account for buoyancy effects. Thermal diffusivities are related to the modelled eddy-viscosity
through the turbulent Prandtl number which takes a typical value of ``Pr_{t}= 1/3`` such that ``\kappa_{2} = \nu/Pr_{t}``.

Tendencies for fourth-order hyperdiffusion are included in the `rhs!` construction, but the
coefficient ``\kappa_{4}`` is ``0`` in this demonstrative case. Hyperdiffusive
tendencies are typically included as a scale-selective diffusion mechanism for high-frequency noise
(e.g. stabilization in GCMs).

Consider components of the viscous stress tensor in three dimensions:

```math
\begin{aligned}
\tau_{xx} = 2\nu \frac{\partial u}{\partial x}, \\

\tau_{yy} = 2\nu \frac{\partial v}{\partial y}, \\

\tau_{zz} = 2\nu \frac{\partial w}{\partial z}, \\

\tau_{xy} = \nu \Big(\frac{\partial u}{\partial y} +  \frac{\partial v}{\partial x}\Big), \\

\tau_{xz} = \nu \Big(\frac{\partial u}{\partial z} +  \frac{\partial w}{\partial x}\Big), \\

\tau_{yz} = \nu \Big(\frac{\partial v}{\partial z} +  \frac{\partial w}{\partial y}\Big).
\end{aligned}
```

Assume terms in the ``y``-direction are neglected (2-dimensional simplicfication). The contributions to the momentum equation are then given by:
```math
\begin{aligned}
(\rho u):  \partial_{x} (\rho \tau_{xx}) + \partial_{z}(\rho\tau_{xz})  &= \partial_x  \Big(2\nu \frac{\partial u}{\partial x}\Big) + \partial_z\Big(\nu \frac{\partial u}{\partial z}\Big) + \partial_z\Big(\nu \frac{\partial w}{\partial x}\Big), \\
(\rho w): \partial_{x} (\rho \tau_{zx})+ \partial_{z}(\rho\tau_{zz})  &= \partial_x\Big(\nu \frac{\partial u}{\partial z}\Big) +  \partial_x\Big(\nu \frac{\partial w}{\partial x}\Big) + \partial_z\Big(2\nu\frac{\partial w}{\partial z} \Big). \\
\end{aligned}
```
Which can be interpreted as, for horizontal-momentum:
1) Horizontal divergence of vertical gradients of cell-centered variables ``u``
2) Vertical divergence of vertical gradients of cell-centered variables ``u``
3) Vertical divergence of horizontal gradients of cell-face variables ``w``

and for vertical-momentum, as:
1) Horizontal divergence of vertical gradients of cell-centered variables ``u``
2) Horizontal divergence of horizontal gradients of cell-face variables ``w``
3) Vertical divergence of vertical gradients of cell-face variables ``w``.

=#

# ## Model Code
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", ".."))

using Test
using StaticArrays, IntervalSets, LinearAlgebra

import ClimaCore: ClimaCore, slab, Spaces, Domains, Meshes, Geometry, Topologies, Spaces, Fields, Operators
using ClimaCore.Geometry
using ClimaCore.Utilities: PlusHalf

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using ClimaCoupler

# Load coupled simulation code
include("../CoupledSims/coupled_sim.jl")

## set up function space
function hvspace_2D(xlim = (-π, π), zlim = (0, 4π), helem = 20, velem = 20, npoly = 1)
    FT = Float64
    vertdomain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_tags = (:bottom, :top),
    )
    vertmesh = Meshes.IntervalMesh(vertdomain, nelems = velem)
    vert_center_space = Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = Domains.IntervalDomain(Geometry.XPoint{FT}(xlim[1]) .. Geometry.XPoint{FT}(xlim[2]), periodic = true)
    horzmesh = Meshes.IntervalMesh(horzdomain; nelems = helem)
    horztopology = Topologies.IntervalTopology(horzmesh)

    quad = Spaces.Quadratures.GLL{npoly + 1}()
    horzspace = Spaces.SpectralElementSpace1D(horztopology, quad)

    hv_center_space = Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)
    return (hv_center_space, hv_face_space)
end

function pressure(ρθ)
    if ρθ >= 0
        return MSLP * (R_d * ρθ / MSLP)^γ
    else
        return NaN
    end
end

Φ(z) = grav * z

abstract type BCtag end
struct ZeroFlux <: BCtag end

bc_divF2C_bottom!(::ZeroFlux, dY, Y, p, t) = Operators.SetValue(Geometry.WVector(0.0))
bc_divF2C_top!(::ZeroFlux, dY, Y, p, t) = Operators.SetValue(Geometry.WVector(0.0))

function init_sea_breeze_2d(x, z)
    θ₀ = atm_T_ini
    cp_d = C_p
    cv_d = C_v
    p₀ = MSLP
    g = grav
    γ = cp_d / cv_d
    z_c = 100.0
    θ_b = atm_T_ini
    θ_p = z < z_c ? rand() - 0.5 : 0.0 # potential temperature perturbation
    θ = θ_b + θ_p # potential temperature
    π_exn = 1.0 - g * z / cp_d / θ # exner function
    T = π_exn * θ # temperature
    p = p₀ * π_exn^(cp_d / R_d) # pressure
    ρ = p / R_d / T # density
    ρθ = ρ * θ # potential temperature density
    return (ρ = ρ, ρθ = ρθ, ρuₕ = ρ * Geometry.UVector(0.0))
end

function atm_rhs!(dY, Y, params, t)
    ρw = Y.ρw
    Yc = Y.Yc
    dYc = dY.Yc
    dρw = dY.ρw

    center_coords = Fields.coordinate_field(axes(Yc))

    ## spectral horizontal operators
    hdiv = Operators.Divergence()
    hgrad = Operators.Gradient()
    hwdiv = Operators.WeakDivergence()
    hwgrad = Operators.WeakGradient()

    ## vertical FD operators with BC's
    vdivf2c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.WVector(0.0)),
        top = Operators.SetValue(Geometry.WVector(0.0)),
    )
    vvdivc2f = Operators.DivergenceC2F(
        bottom = Operators.SetDivergence(Geometry.WVector(0.0)),
        top = Operators.SetDivergence(Geometry.WVector(0.0)),
    )
    uvdivf2c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.WVector(0.0) ⊗ Geometry.UVector(0.0)),
        top = Operators.SetValue(Geometry.WVector(0.0) ⊗ Geometry.UVector(0.0)),
    )
    If = Operators.InterpolateC2F(bottom = Operators.Extrapolate(), top = Operators.Extrapolate())
    Ic = Operators.InterpolateF2C()
    ∂ = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.WVector(0.0)),
        top = Operators.SetValue(Geometry.WVector(0.0)),
    )
    ∂f = Operators.GradientC2F()
    ∂c = Operators.GradientF2C()
    B = Operators.SetBoundaryOperator(
        bottom = Operators.SetValue(Geometry.WVector(0.0)),
        top = Operators.SetValue(Geometry.WVector(0.0)),
    )

    ∇_z_ρθ = Operators.DivergenceF2C(
        bottom = bc_divF2C_bottom!(params.bc.ρθ.bottom, dY, Y, params, t),
        top = bc_divF2C_top!(params.bc.ρθ.top, dY, Y, params, t),
    )

    uₕ = @. Yc.ρuₕ / Yc.ρ
    w = @. ρw / If(Yc.ρ)
    wc = @. Ic(ρw) / Yc.ρ
    p = @. pressure(Yc.ρθ)
    θ = @. Yc.ρθ / Yc.ρ
    Yfρ = @. If(Yc.ρ)

    ### HYPERVISCOSITY
    ## 1) compute hyperviscosity coefficients
    @. dYc.ρθ = hwdiv(hgrad(θ))
    @. dYc.ρuₕ = hwdiv(hgrad(uₕ))
    @. dρw = hwdiv(hgrad(w))
    Spaces.weighted_dss!(dYc)
    Spaces.weighted_dss!(dρw)

    κ₄ = 0.0 # m^4/s
    @. dYc.ρθ = -κ₄ * hwdiv(Yc.ρ * hgrad(dYc.ρθ))
    @. dYc.ρuₕ = -κ₄ * hwdiv(Yc.ρ * hgrad(dYc.ρuₕ))
    @. dρw = -κ₄ * hwdiv(Yfρ * hgrad(dρw))

    ## density
    @. dYc.ρ = -∂(ρw)
    @. dYc.ρ -= hdiv(Yc.ρuₕ)

    ## potential temperature
    @. dYc.ρθ += -(∇_z_ρθ(ρw * If(Yc.ρθ / Yc.ρ)))
    @. dYc.ρθ -= hdiv(uₕ * Yc.ρθ)

    ## horizontal momentum
    Ih = Ref(Geometry.Axis2Tensor((Geometry.UAxis(), Geometry.UAxis()), @SMatrix [1.0]))
    @. dYc.ρuₕ += -uvdivf2c(ρw ⊗ If(uₕ))
    @. dYc.ρuₕ -= hdiv(Yc.ρuₕ ⊗ uₕ + p * Ih)

    ## vertical momentum
    @. dρw +=
        B(Geometry.transform(Geometry.WAxis(), -(∂f(p)) - If(Yc.ρ) * ∂f(Φ(center_coords.z))) - vvdivc2f(Ic(ρw ⊗ w)))
    uₕf = @. If(Yc.ρuₕ / Yc.ρ) # requires boundary conditions
    @. dρw -= hdiv(uₕf ⊗ ρw)

    ## DIFFUSION
    κ₂ = 5.0 # m^2/s
    ##  1a) horizontal div of horizontal grad of horiz momentun
    @. dYc.ρuₕ += hwdiv(κ₂ * (Yc.ρ * hgrad(Yc.ρuₕ / Yc.ρ)))
    ##  1b) vertical div of vertical grad of horiz momentun
    @. dYc.ρuₕ += uvdivf2c(κ₂ * (Yfρ * ∂f(Yc.ρuₕ / Yc.ρ)))

    ##  1c) horizontal div of horizontal grad of vert momentum
    @. dρw += hwdiv(κ₂ * (Yfρ * hgrad(ρw / Yfρ)))
    ##  1d) vertical div of vertical grad of vert momentun
    @. dρw += vvdivc2f(κ₂ * (Yc.ρ * ∂c(ρw / Yfρ)))

    ##  2a) horizontal div of horizontal grad of potential temperature
    @. dYc.ρθ += hwdiv(κ₂ * (Yc.ρ * hgrad(Yc.ρθ / Yc.ρ)))
    ##  2b) vertical div of vertial grad of potential temperature
    @. dYc.ρθ += ∇_z_ρθ(κ₂ * (Yfρ * ∂f(Yc.ρθ / Yc.ρ)))

    Spaces.weighted_dss!(dYc)
    Spaces.weighted_dss!(dρw)
    return dY
end

## init simulation
function atm_init(; xmin = -500, xmax = 500, zmin = 0, zmax = 1000, npoly = 3, helem = 20, velem = 20, bc = nothing)

    ## construct domain spaces
    hv_center_space, hv_face_space = hvspace_2D((xmin, xmax), (zmin, zmax), helem, velem, npoly) # [m]
    center_coords = Fields.coordinate_field(hv_center_space)
    face_coords = Fields.coordinate_field(hv_face_space)
    domain = (hv_center_space = hv_center_space, hv_face_space = hv_face_space)

    ## initialize prognostic variables
    Yc = map(center_coords) do coord
        sea_breeze = init_sea_breeze_2d(coord.x, coord.z)
        sea_breeze
    end

    ρw = map(face_coords) do coord
        Geometry.WVector(0.0)
    end

    Y = Fields.FieldVector(Yc = Yc, ρw = ρw)

    ## select boundary conditions
    if bc === nothing
        bc = (
            ρθ = (bottom = CoupledFlux(), top = ZeroFlux()),
            ρu = nothing, # for now BCs are hard coded, except for ρθ
        )
    end

    return Y, bc, domain
end

# ## Coupled Atmos Wrappers
## Atmos Simulation - later to live in ClimaAtmos
struct AtmosSim <: AbstractAtmosSim
    integrator::Any
end

function AtmosSim(Y_init, t_start, dt, t_end, timestepper, p, saveat, callbacks = CallbackSet())
    ode_algo = CTS.ExplicitAlgorithm(timestepper)
    ode_function = CTS.ClimaODEFunction(T_exp! = atm_rhs!)

    problem = ODEProblem(ode_function, Y_init, (t_start, t_end), p)
    atm_integ = init(
        problem,
        ode_algo,
        dt = dt,
        saveat = saveat,
        adaptive = false,
        progress = true,
        progress_message = (dt, u, params, t) -> t,
        callback = callbacks,
    )

    return AtmosSim(atm_integ)
end

function coupler_push!(coupler::CouplerState, atmos::AtmosSim)
    coupler_put!(coupler, :F_sfc, atmos.integrator.u.F_sfc, atmos)
end

function coupler_pull!(atmos::AtmosSim, coupler::CouplerState)
    ## reset flux accumulator
    atmos.integrator.u.F_sfc .= 0.0 # reset surface flux to be accumulated

    T_sfc_ocean = coupler_get(coupler, :T_sfc_ocean, atmos)
    T_sfc_land = coupler_get(coupler, :T_sfc_land, atmos)
    atmos.integrator.p.T_sfc .= T_sfc_land .+ T_sfc_ocean
end

#=
## Coupled Boundary Conditions

The standalone atmosphere model uses two boundary condition methods in its tendency:
`bc_divF2C_bottom!` and `bc_divF2C_top!`. Since the bottom boundary is coupled, `bc_divF2C_bottom!`
must be altered when running in coupled mode to properly calculate and accumulate the boundary flux
from the ocean and land components.

To solve this, a `CoupledFlux` boundary tag is set for the bottom boundary during initialization.
Then, a new method of `bc_divF2C_bottom!` is written to dispatch on the `CoupledFlux` boundary tag.
This method can then compute the flux appropriately.
=#
struct CoupledFlux <: BCtag end
function bc_divF2C_bottom!(::CoupledFlux, dY, Y, p, t)
    ## flux calculation
    Yc = Y.Yc
    uₕ = Yc.ρuₕ ./ Yc.ρ
    ρw = Y.ρw
    If2c = Operators.InterpolateF2C()
    Ic2f = Operators.InterpolateC2F(bottom = Operators.Extrapolate(), top = Operators.Extrapolate())
    w = If2c.(ρw) ./ Yc.ρ
    cuv = @. Geometry.UWVector(uₕ)
    windspeed = @. norm(cuv)
    windspeed_boundary = Fields.level(windspeed, 1)
    θ_boundary = Fields.level(Yc.ρθ ./ Yc.ρ, 1)
    ρ_boundary = Fields.level(Yc.ρ, 1)

    ## build atmos face fields on surface boundary space to enable broadcasting
    windspeed_boundary = Fields.Field(Fields.field_values(windspeed_boundary), axes(p.T_sfc))
    θ_boundary = Fields.Field(Fields.field_values(θ_boundary), axes(p.T_sfc))
    ρ_boundary = Fields.Field(Fields.field_values(ρ_boundary), axes(p.T_sfc))

    λ = @. p.cpl_p.C_p * p.cpl_p.C_H * ρ_boundary * windspeed_boundary
    dθ = @. θ_boundary - p.T_sfc
    heat_flux = @. -λ * dθ
    @. dY.F_sfc += heat_flux # accumulation

    return Operators.SetValue(Geometry.WVector.(heat_flux))
end
