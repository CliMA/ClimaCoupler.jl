include("boundary_conditions_adapted.jl")
"""
    calc_boundary_state!(::NumericalFluxSecondOrder, ::Impenetrable{FreeSlip}, ::ModelSetup)
apply free slip boundary condition for velocity
apply zero numerical flux in the normal direction
"""
function calc_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Impenetrable{FreeSlip},
    ::Union{ModelSetup},
    ::Nothing,
    state⁺,
    gradflux⁺,
    hyperflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    hyperflux⁻,
    aux⁻,
    t,
    _...,
)
    state⁺.ρu = state⁻.ρu

    return nothing
end

function calc_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Impenetrable{FreeSlip},
    ::Union{ModelSetup},
    ::ConstantViscosity,
    state⁺,
    gradflux⁺,
    hyperflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    hyperflux⁻,
    aux⁻,
    t,
    _...,
)
    state⁺.ρu = state⁻.ρu
    gradflux⁺.ν∇u = n⁻ * (@SVector [-0, -0, -0])'

    return nothing
end

"""
    calc_boundary_state!(::NumericalFluxSecondOrder, ::Impenetrable{NoSlip}, ::ModelSetup)
apply no slip boundary condition for velocity
sets ghost point to have no numerical flux on the boundary for U
"""
@inline function calc_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Impenetrable{NoSlip},
    ::Union{ModelSetup},
    ::Nothing,
    state⁺,
    gradflux⁺,
    hyperflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    hyperflux⁻,
    aux⁻,
    t,
    _...,
)
    state⁺.ρu = -state⁻.ρu

    return nothing
end

@inline function calc_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Impenetrable{NoSlip},
    ::Union{ModelSetup},
    ::ConstantViscosity,
    state⁺,
    gradflux⁺,
    hyperflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    hyperflux⁻,
    aux⁻,
    t,
    _...,
)
    state⁺.ρu = -state⁻.ρu
    gradflux⁺.ν∇u = gradflux⁻.ν∇u

    return nothing
end

"""
    calc_boundary_state!(::NumericalFluxSecondOrder, ::Insulating, ::HBModel)

apply insulating boundary condition for velocity
sets ghost point to have no numerical flux on the boundary for κ∇θ
"""
@inline function calc_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Union{Insulating, CoupledSecondaryBoundary, CoupledPrimaryBoundary},
    ::Union{ModelSetup},
    ::Nothing,
    state⁺,
    gradflux⁺,
    hyperflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    hyperflux⁻,
    aux⁻,
    t,
    _...,
)
    state⁺.ρθ = state⁻.ρθ

    return nothing
end

@inline function calc_boundary_state!(
    ::NumericalFluxSecondOrder,
    ::Union{Insulating, CoupledSecondaryBoundary, CoupledPrimaryBoundary},
    ::Union{ModelSetup},
    ::ConstantViscosity,
    state⁺,
    gradflux⁺,
    hyperflux⁺,
    aux⁺,
    n⁻,
    state⁻,
    gradflux⁻,
    hyperflux⁻,
    aux⁻,
    t,
    _...,
)
    state⁺.ρθ = state⁻.ρθ
    gradflux⁺.κ∇θ = n⁻ * -0

    return nothing
end
